# Redaction note (2026-09-25)

`queries.log` and `raw/provenance.jsonl` contained two MAST download URLs that redirected to Amazon S3 with
time-limited *presigned* authorisation parameters (`X-Amz-Algorithm`, `X-Amz-Date`, `X-Amz-SignedHeaders`,
`X-Amz-Credential`, `X-Amz-Expires`, `X-Amz-Signature`). These parameters belonged to the MAST service, not to the
author, and expired one hour after issue. Their values were replaced by `REDACTED-presigned-auth-2026-09-25` before
deposit. Nothing else in the two files was changed: the MAST product identifiers, the request times, the check results
and the recorded checksums are intact. The base URLs (without the authorisation parameters) still identify the products.
