# GCP Cloud Run Deployment Guide: ExoElenchus Portal

This guide explains how to deploy the **TIC 4206066 Discovery & Verification Portal** to Google Cloud Platform (GCP) Cloud Run.

## 1. Prerequisites
- Google Cloud SDK (`gcloud` CLI installed and authenticated)
- A GCP project with Billing enabled
- Cloud Run API and Cloud Build API enabled:
  ```bash
  gcloud services enable run.googleapis.com cloudbuild.googleapis.com
  ```

---

## 2. Quick Deploy via Script

Run the automated deployment script from the repository root:
```bash
./deploy_cloud_run.sh
```

Or configure custom parameters:
```bash
SERVICE_NAME="tic4206066-portal" REGION="us-central1" PROJECT_ID="my-gcp-project" ./deploy_cloud_run.sh
```

---

## 3. Direct gcloud Command

Deploy directly from source with zero build setup:
```bash
gcloud run deploy tic4206066-portal \
  --source . \
  --region us-central1 \
  --platform managed \
  --allow-unauthenticated \
  --memory 1Gi \
  --cpu 1 \
  --port 8080
```

---

## 4. Automated CI/CD via Cloud Build

Submit the build trigger using `cloudbuild.yaml`:
```bash
gcloud builds submit --config=cloudbuild.yaml
```

---

## 5. Local Docker Testing

Build and test the container locally before deploying to GCP:
```bash
# Build Docker image
docker build -t tic4206066-portal .

# Run container locally on port 8080
docker run -p 8080:8080 -e PORT=8080 tic4206066-portal

# Access dashboard
open http://localhost:8080
```

---

## 6. Production Architecture
- **Web Frontend:** Stitch-designed astrophysics command deck with live countdown ticker and real-time telemetry console.
- **Backend API:** FastAPI ASGI service providing `/health`, `/api/astrophysics/tic4206066/scheduler`, `/api/astrophysics/tic4206066/verify`, `/api/astrophysics/tic4206066/status`, `/api/astrophysics/tic4206066/reproduce`.
- **Runtime:** Fully containerized on Python 3.12-slim with stateless auto-scaling (0 to N instances) and sub-second cold starts.
