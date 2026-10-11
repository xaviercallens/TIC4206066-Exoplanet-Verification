"""
Test Suite: FastAPI Backend & GCP Cloud Run Portal for TIC 4206066.
Verifies all REST API endpoints:
- GET /: Dashboard HTML delivery
- GET /health: Health check
- GET /api/astrophysics/tic4206066/scheduler: Live countdown and scheduler status
- GET /api/astrophysics/tic4206066/verify: Full verification workflow output and Lean 4 proof status
- POST /api/astrophysics/tic4206066/scheduler/dry-run: Dry-run workflow execution
- GET /api/astrophysics/tic4206066/status: Certified receipt & target metadata
"""

import pytest
from fastapi.testclient import TestClient

from app import app

client = TestClient(app)


def test_health_endpoint():
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "healthy"
    assert "service" in data


def test_dashboard_delivery():
    response = client.get("/")
    assert response.status_code == 200
    assert "text/html" in response.headers["content-type"]
    text = response.text
    assert "TIC 4206066" in text
    assert "ExoElenchus" in text


def test_scheduler_status():
    response = client.get("/api/astrophysics/tic4206066/scheduler")
    assert response.status_code == 200
    data = response.json()
    assert "time_counter" in data
    tc = data["time_counter"]
    assert "target_date_utc" in tc
    assert "current_time_utc" in tc
    assert "formatted_countdown" in tc
    assert "days" in tc


def test_verification_endpoint():
    response = client.get("/api/astrophysics/tic4206066/verify")
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert "provenance_hash" in data
    assert "lean4_proof" in data
    assert data["lean4_proof"]["status"] in ("ELABORATED_ZERO_SORRY", "UNKNOWN")
    assert "workflow_summary" in data
    steps = data["workflow_summary"]["steps"]
    assert len(steps) >= 7


def test_dry_run_endpoint():
    response = client.post("/api/astrophysics/tic4206066/scheduler/dry-run")
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["is_dry_run"] is True


def test_status_endpoint():
    response = client.get("/api/astrophysics/tic4206066/status")
    assert response.status_code == 200
    data = response.json()
    assert "target_id" in data or "target" in data


def test_lean_endpoint():
    response = client.get("/api/astrophysics/tic4206066/lean")
    assert response.status_code == 200
    data = response.json()
    assert "status" in data
    assert "theorems" in data["status"]
    assert len(data["status"]["theorems"]) >= 6
    assert data["status"]["zero_sorry"] is True


def test_reproduce_endpoint():
    response = client.post("/api/astrophysics/tic4206066/reproduce")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"
    assert data["target"] == "TIC 4206066"
    assert "data_integrity" in data

