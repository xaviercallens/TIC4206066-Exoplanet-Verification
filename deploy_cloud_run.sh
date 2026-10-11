#!/usr/bin/env bash
# Deploy ExoElenchus TIC 4206066 Discovery Portal to Google Cloud Run
# Optimized for low cost (scale to zero, cpu-throttling, memory 512Mi, max-instances 2)
set -euo pipefail

SERVICE_NAME="${SERVICE_NAME:-exoplanet-tic4206066-portal}"
REGION="${REGION:-europe-west1}"
PROJECT_ID="${PROJECT_ID:-$(gcloud config get-value project 2>/dev/null || echo 'gen-lang-client-0625573011')}"

if [[ -z "$PROJECT_ID" ]]; then
  echo "Error: No active GCP project configured. Set PROJECT_ID or run: gcloud config set project <PROJECT_ID>" >&2
  exit 1
fi

echo "============================================================"
echo " Deploying TIC 4206066 Discovery Portal to Google Cloud Run"
echo " Service:        $SERVICE_NAME"
echo " Region:         $REGION"
echo " Project:        $PROJECT_ID"
echo " Cost Profile:   Scale-to-zero (min: 0, max: 2), CPU throttled"
echo "============================================================"

gcloud run deploy "$SERVICE_NAME" \
  --source . \
  --project "$PROJECT_ID" \
  --region "$REGION" \
  --platform managed \
  --allow-unauthenticated \
  --memory 512Mi \
  --cpu 1 \
  --min-instances 0 \
  --max-instances 2 \
  --concurrency 80 \
  --cpu-throttling \
  --timeout 60s \
  --port 8080

echo ""
echo "Deployment successful! Service URL:"
gcloud run services describe "$SERVICE_NAME" --project "$PROJECT_ID" --region "$REGION" --format 'value(status.url)'
