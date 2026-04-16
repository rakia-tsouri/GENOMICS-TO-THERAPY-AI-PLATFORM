import httpx
from typing import Optional
from app.utils.cache import cache, TTL_7D
from app.models import ProteinAnnotation, ActiveSite, BindingSite
import logging

logger = logging.getLogger(__name__)

async def get_uniprot_annotation(uniprot_id: str) -> Optional[ProteinAnnotation]:
    """
    Fetch protein annotation from UniProt REST API.
    """
    if not uniprot_id:
        return None

    # Check cache
    cache_key = f"uniprot:{uniprot_id}"
    cached_result = cache.get(cache_key)
    if cached_result:
        return ProteinAnnotation(**cached_result)

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            response = await client.get(url)
            if response.status_code != 200:
                logger.warning(f"UniProt API returned {response.status_code} for {uniprot_id}")
                return None
            
            data = response.json()
            
            # Extract basic function
            function = None
            if "comments" in data:
                for comment in data["comments"]:
                    if comment.get("commentType") == "FUNCTION":
                        function = comment.get("texts", [{}])[0].get("value")
                        break
            
            # Extract domains, active sites, binding sites
            domains = []
            active_sites = []
            binding_sites = []
            if "features" in data:
                for feature in data["features"]:
                    f_type = feature.get("type")
                    description = feature.get("description", "")
                    location = feature.get("location", {})
                    
                    if f_type == "Domain":
                        domains.append(description)
                    elif f_type == "Active site":
                        start = location.get("start", {}).get("value")
                        if start:
                            active_sites.append(ActiveSite(position=start, description=description))
                    elif f_type == "Binding site":
                        start = location.get("start", {}).get("value")
                        end = location.get("end", {}).get("value")
                        ligand = feature.get("ligand", {}).get("name", "Unknown")
                        if start and end:
                            binding_sites.append(BindingSite(positions=[start, end], ligand=ligand))
                        elif start:
                            binding_sites.append(BindingSite(positions=[start], ligand=ligand))

            # Extract diseases
            diseases = []
            if "comments" in data:
                for comment in data["comments"]:
                    if comment.get("commentType") == "DISEASE":
                        disease_obj = comment.get("disease", {})
                        if "diseaseId" in disease_obj:
                            diseases.append(disease_obj.get("diseaseId"))

            result = ProteinAnnotation(
                function=function,
                domains=domains,
                active_sites=active_sites,
                binding_sites=binding_sites,
                diseases=diseases
            )

            # Save to cache
            cache.set(cache_key, result.model_dump(), ttl=TTL_7D)
            return result

    except Exception as e:
        logger.error(f"UniProt API request failed for {uniprot_id}: {e}")
        return None
