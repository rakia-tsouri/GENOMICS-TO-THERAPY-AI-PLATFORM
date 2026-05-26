from fastapi import APIRouter, Depends, HTTPException, status
from fastapi.security import OAuth2PasswordRequestForm
from sqlalchemy.orm import Session

from ..database import get_db
from ..deps import get_current_user
from ..models import User
from ..schemas.auth import LoginRequest, RegisterRequest, Token
from ..schemas.user import UserRead
from ..security import create_access_token, hash_password, verify_password

router = APIRouter(prefix="/auth", tags=["auth"])


def _authenticate(db: Session, email: str, password: str) -> User:
    user = db.query(User).filter(User.email == email).first()
    if not user or not verify_password(password, user.hashed_password):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="Incorrect email or password"
        )
    return user


@router.post("/register", response_model=Token, status_code=201)
def register(body: RegisterRequest, db: Session = Depends(get_db)):
    if db.query(User).filter(User.email == body.email).first():
        raise HTTPException(status_code=400, detail="Email already registered")
    # First ever user becomes admin; everyone else is a researcher.
    role = "admin" if db.query(User).count() == 0 else "researcher"
    user = User(
        email=body.email,
        full_name=body.full_name,
        hashed_password=hash_password(body.password),
        role=role,
    )
    db.add(user)
    db.commit()
    db.refresh(user)
    return Token(access_token=create_access_token(str(user.id), {"role": user.role}))


@router.post("/login", response_model=Token)
def login(body: LoginRequest, db: Session = Depends(get_db)):
    user = _authenticate(db, body.email, body.password)
    return Token(access_token=create_access_token(str(user.id), {"role": user.role}))


@router.post("/token", response_model=Token)
def token(form: OAuth2PasswordRequestForm = Depends(), db: Session = Depends(get_db)):
    """OAuth2 password flow (used by Swagger 'Authorize' and as a generic token endpoint)."""
    user = _authenticate(db, form.username, form.password)
    return Token(access_token=create_access_token(str(user.id), {"role": user.role}))


@router.get("/me", response_model=UserRead)
def me(user: User = Depends(get_current_user)):
    return user
