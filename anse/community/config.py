"""
ANSE Community Engine - Configuration and Target Domains
"""

import os
from dataclasses import dataclass, field
from enum import StrEnum
from typing import Any


class ResearchDomain(StrEnum):
    RUST_LINUX = "rust_linux"
    AI_OPTIMIZATION = "ai_optimization"
    OPEN_WEIGHTS = "open_weights"
    ASTROPHYSICS = "astrophysics"
    QUANTUM_COMPUTING = "quantum_computing"
    QUANTUM_PHYSICS = "quantum_physics"
    HPC = "hpc"


DOMAIN_SUBREDDITS: dict[ResearchDomain, list[str]] = {
    ResearchDomain.RUST_LINUX: ["rust", "linux", "kernel", "systems"],
    ResearchDomain.AI_OPTIMIZATION: ["MachineLearning", "LocalLLaMA", "deeplearning", "ArtificialInteligence"],
    ResearchDomain.OPEN_WEIGHTS: ["LocalLLaMA", "OpenSourceAI", "MachineLearning"],
    ResearchDomain.ASTROPHYSICS: ["Astrophysics", "space", "Cosmology", "astronomy"],
    ResearchDomain.QUANTUM_COMPUTING: ["QuantumComputing", "Quantum"],
    ResearchDomain.QUANTUM_PHYSICS: ["Physics", "TheoreticalPhysics", "Quantum"],
    ResearchDomain.HPC: ["HPC", "supercomputing", "parallelprogramming"],
}

DOMAIN_KEYWORDS: dict[ResearchDomain, list[str]] = {
    ResearchDomain.RUST_LINUX: [
        "rust", "linux", "kernel", "ebpf", "io_uring", "memory safety",
        "simd", "zero-cost", "concurrency", "allocator"
    ],
    ResearchDomain.AI_OPTIMIZATION: [
        "quantization", "gguf", "awq", "flashattention", "lora", "unsloth",
        "parameter budget", "jepa", "energy minimization", "zero allocation",
        "kernel optimization", "inference latency"
    ],
    ResearchDomain.OPEN_WEIGHTS: [
        "open weights", "open source ai", "reproducible", "local llm",
        "ollama", "vllm", "open-r1", "checkpoint", "leaderboard"
    ],
    ResearchDomain.ASTROPHYSICS: [
        "desi", "bao", "baryon acoustic oscillations", "eboss", "cosmology",
        "flcdm", "h0 tension", "dark energy", "gravitational waves"
    ],
    ResearchDomain.QUANTUM_COMPUTING: [
        "quantum computing", "qubit", "tqec", "surface code", "quantum error correction",
        "clifford", "tensor network", "quantum circuit"
    ],
    ResearchDomain.QUANTUM_PHYSICS: [
        "symplectic", "hamiltonian", "lattice gauge", "quantum field theory",
        "conservation law", "lie algebra", "unitarity"
    ],
    ResearchDomain.HPC: [
        "hpc", "mpi", "avx-512", "roofline", "cuda", "distributed training",
        "numa", "cache locality", "bandwidth"
    ],
}

BANNED_BUZZWORDS: list[str] = [
    "game-changing", "revolutionary", "groundbreaking", "paradigm shift",
    "unprecedented", "mind-blowing", "miraculous", "disruptive",
    "quantum leap", "holy grail"
]


@dataclass
class CommunityConfig:
    # Reddit API Credentials
    reddit_client_id: str = field(default_factory=lambda: os.getenv("REDDIT_CLIENT_ID", ""))
    reddit_client_secret: str = field(default_factory=lambda: os.getenv("REDDIT_CLIENT_SECRET", ""))
    reddit_user_agent: str = field(
        default_factory=lambda: os.getenv("REDDIT_USER_AGENT", "ANSE-Community-Agent/1.0 (academic research dissemination)")
    )
    reddit_username: str = field(default_factory=lambda: os.getenv("REDDIT_USERNAME", ""))
    reddit_password: str = field(default_factory=lambda: os.getenv("REDDIT_PASSWORD", ""))
    reddit_base_url: str = "https://www.reddit.com"

    # API Endpoints
    zenodo_api_base: str = "https://zenodo.org/api/records"
    arxiv_api_base: str = "http://export.arxiv.org/api/query"

    # Ethical & Community Guardrails
    max_promotional_ratio: float = 0.10  # 9:1 community ratio (Reddit self-promotion guideline)
    require_submission_statement: bool = True
    require_human_approval: bool = True
    min_confidence_score: float = 0.92

    # Local Repository Roots
    papers_dir: str = "papers"
    results_dir: str = "results"
    formal_dir: str = "formal"

    def to_dict(self) -> dict[str, Any]:
        return {
            "reddit_configured": bool(self.reddit_client_id and self.reddit_client_secret),
            "max_promotional_ratio": self.max_promotional_ratio,
            "require_submission_statement": self.require_submission_statement,
            "require_human_approval": self.require_human_approval,
            "papers_dir": self.papers_dir,
        }
