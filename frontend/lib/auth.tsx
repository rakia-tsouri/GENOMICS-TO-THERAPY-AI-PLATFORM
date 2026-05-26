"use client";

import {
  createContext,
  useCallback,
  useContext,
  useEffect,
  useMemo,
  useState,
  type ReactNode,
} from "react";
import {
  clearToken,
  get,
  post,
  setToken,
  getToken,
} from "./api";
import type { AuthToken, User } from "./types";

interface AuthContextValue {
  user: User | null;
  loading: boolean;
  isAuthenticated: boolean;
  isAdmin: boolean;
  login: (email: string, password: string) => Promise<User>;
  register: (
    email: string,
    password: string,
    fullName: string
  ) => Promise<User>;
  logout: () => void;
  refresh: () => Promise<void>;
  setUser: (user: User) => void;
}

const AuthContext = createContext<AuthContextValue | undefined>(undefined);

export function AuthProvider({ children }: { children: ReactNode }) {
  const [user, setUserState] = useState<User | null>(null);
  const [loading, setLoading] = useState(true);

  const loadMe = useCallback(async () => {
    const token = getToken();
    if (!token) {
      setUserState(null);
      setLoading(false);
      return;
    }
    try {
      const me = await get<User>("/auth/me");
      setUserState(me);
    } catch {
      clearToken();
      setUserState(null);
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    void loadMe();
  }, [loadMe]);

  const login = useCallback(
    async (email: string, password: string) => {
      const res = await post<AuthToken>("/auth/login", { email, password });
      setToken(res.access_token);
      const me = await get<User>("/auth/me");
      setUserState(me);
      setLoading(false);
      return me;
    },
    []
  );

  const register = useCallback(
    async (email: string, password: string, fullName: string) => {
      const res = await post<AuthToken>("/auth/register", {
        email,
        password,
        full_name: fullName,
      });
      setToken(res.access_token);
      const me = await get<User>("/auth/me");
      setUserState(me);
      setLoading(false);
      return me;
    },
    []
  );

  const logout = useCallback(() => {
    clearToken();
    setUserState(null);
  }, []);

  const refresh = useCallback(async () => {
    await loadMe();
  }, [loadMe]);

  const setUser = useCallback((u: User) => {
    setUserState(u);
  }, []);

  const value = useMemo<AuthContextValue>(
    () => ({
      user,
      loading,
      isAuthenticated: !!user,
      isAdmin: user?.role === "admin",
      login,
      register,
      logout,
      refresh,
      setUser,
    }),
    [user, loading, login, register, logout, refresh, setUser]
  );

  return <AuthContext.Provider value={value}>{children}</AuthContext.Provider>;
}

export function useAuth(): AuthContextValue {
  const ctx = useContext(AuthContext);
  if (!ctx) {
    throw new Error("useAuth must be used within an AuthProvider");
  }
  return ctx;
}
