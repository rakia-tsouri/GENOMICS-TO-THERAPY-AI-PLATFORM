import os
import zlib
import urllib.request

def encode_6bit(b):
    if b < 10:
        return chr(48 + b)
    b -= 10
    if b < 26:
        return chr(65 + b)
    b -= 26
    if b < 26:
        return chr(97 + b)
    b -= 26
    if b == 0:
        return '-'
    if b == 1:
        return '_'
    return '?'

def encode_3bytes(b1, b2, b3):
    c1 = b1 >> 2
    c2 = ((b1 & 0x3) << 4) | (b2 >> 4)
    c3 = ((b2 & 0xF) << 2) | (b3 >> 6)
    c4 = b3 & 0x3F
    return encode_6bit(c1) + encode_6bit(c2) + encode_6bit(c3) + encode_6bit(c4)

def encode_plantuml(text):
    zlib_comp = zlib.compressobj(9, zlib.DEFLATED, -zlib.MAX_WBITS, zlib.DEF_MEM_LEVEL, 0)
    compressed = zlib_comp.compress(text.encode('utf-8')) + zlib_comp.flush()
    res = ""
    i = 0
    while i < len(compressed):
        if i + 2 < len(compressed):
            res += encode_3bytes(compressed[i], compressed[i+1], compressed[i+2])
        elif i + 1 < len(compressed):
            b1, b2 = compressed[i], compressed[i+1]
            c1 = b1 >> 2
            c2 = ((b1 & 0x3) << 4) | (b2 >> 4)
            c3 = (b2 & 0xF) << 2
            res += encode_6bit(c1) + encode_6bit(c2) + encode_6bit(c3)
        else:
            b1 = compressed[i]
            c1 = b1 >> 2
            c2 = (b1 & 0x3) << 4
            res += encode_6bit(c1) + encode_6bit(c2)
        i += 3
    return res

def main():
    diagrams_dir = r"c:\Users\MSI\GENOMICS-TO-THERAPY-AI-PLATFORM\docs\report\diagrams"
    images_dir = r"c:\Users\MSI\GENOMICS-TO-THERAPY-AI-PLATFORM\docs\report\images"
    os.makedirs(images_dir, exist_ok=True)
    
    files_to_generate = ["activity_genomics.puml"]
    
    for filename in files_to_generate:
        puml_path = os.path.join(diagrams_dir, filename)
        if not os.path.exists(puml_path):
            print(f"File not found: {puml_path}")
            continue
            
        print(f"Processing {filename}...")
        with open(puml_path, "r", encoding="utf-8") as f:
            puml_text = f.read()
        
        encoded = encode_plantuml(puml_text)
        url = f"http://www.plantuml.com/plantuml/png/{encoded}"
        
        try:
            req = urllib.request.Request(url, headers={
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
            })
            with urllib.request.urlopen(req, timeout=30) as response:
                data = response.read()
                png_name = filename.replace(".puml", ".png")
                png_path = os.path.join(images_dir, png_name)
                with open(png_path, "wb") as png_f:
                    png_f.write(data)
                print(f"Saved {png_name} ({len(data)} bytes)")
        except Exception as e:
            print(f"Error fetching {filename}: {e}")
            print(f"URL: {url}")

if __name__ == "__main__":
    main()
