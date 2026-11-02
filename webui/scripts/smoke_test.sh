#!/usr/bin/env bash
# Smoke test of a running web UI stack (`docker compose --env-file oktoberfest-webui.env up -d` in webui/): API, worker,
# frontend, the quantification / localization options and the tools inside the image. Run from webui/:
#   scripts/smoke_test.sh [base url]
# Uses oktoberfest-webui.env for docker compose and API_PORT; override with COMPOSE="docker compose ...".
set -euo pipefail

ENV_FILE=oktoberfest-webui.env
if [[ -z ${API_PORT:-} && -f $ENV_FILE ]]; then
    API_PORT=$(sed -n 's/^API_PORT=//p' "$ENV_FILE" | tail -1)
fi
BASE=${1:-http://localhost:${API_PORT:-8080}}
COMPOSE=${COMPOSE:-docker compose --env-file $ENV_FILE}
TIMEOUT=${SMOKE_TIMEOUT:-180}
COOKIE="oktoberfest_session=smoke-test"  # hosted mode identifies users by this cookie

fail() { echo "FAIL: $*" >&2; exit 1; }
api() { curl -fsS -b "$COOKIE" "$@"; }

echo "== waiting for the API and the worker at $BASE"
for ((i = 0; i < TIMEOUT; i += 5)); do
    health=$(curl -fsS "$BASE/api/v1/health" 2>/dev/null || true)
    [[ $health == *'"worker":"up"'* ]] && break
    sleep 5
done
[[ $health == *'"worker":"up"'* ]] || fail "API/worker not up after ${TIMEOUT}s: ${health:-no response}"
echo "$health"

echo "== frontend"
api "$BASE/" | grep -q 'id="app"' || fail "index.html is not served"

echo "== meta endpoints"
quant=$(api "$BASE/api/v1/meta/quantification")
[[ $quant == *'"flashlfq_available":true'* ]] || fail "FlashLFQ not available: $quant"
api "$BASE/api/v1/meta/localization-ptms" | grep -q '"21"' || fail "phospho missing from localization PTMs"

echo "== job creation with quantification and localization"
create() {
    curl -sS -o /dev/null -w '%{http_code}' -b "$COOKIE" -H 'Content-Type: application/json' \
        -d "$1" "$BASE/api/v1/jobs"
}
job='{"job_type": "Rescoring", "config": {
    "inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"},
    "quantification": true, "fastaDigestOptions": {}}}'
[[ $(create "$job") == 201 ]] || fail "search engine quantification job rejected"
job='{"job_type": "Rescoring", "config": {
    "inputs": {"search_results_type": "Msfragger", "spectra_type": "mzml"}, "fastaDigestOptions": {},
    "quantification": true, "quantificationOptions": {"precursor": "flashlfq"}}}'
[[ $(create "$job") == 201 ]] || fail "flashlfq quantification job rejected"
job='{"job_type": "Rescoring", "config": {
    "inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"},
    "models": {"intensity": "Prosit_2024_intensity_PTMs_gl", "irt": "Prosit_2024_irt_PTMs_gl"},
    "ptm_localization": true, "massTolerance": "auto"}}'
[[ $(create "$job") == 201 ]] || fail "localization job rejected"
job='{"job_type": "Rescoring", "config": {
    "inputs": {"search_results_type": "Maxquant", "spectra_type": "mzml"}, "ptm_localization": true}}'
[[ $(create "$job") == 422 ]] || fail "localization without a PTM model was not rejected"

echo "== tools in the worker image"
$COMPOSE exec -T worker sh -c '
    set -e
    python -c "import oktoberfest, spectrum_fundamentals, spectrum_io, directlfq; print(\"oktoberfest\", oktoberfest.__version__)"
    python -c "from oktoberfest.utils.config import QUANT_PRECURSOR_METHODS; assert \"flashlfq\" in QUANT_PRECURSOR_METHODS"
    percolator --help 2>&1 | grep -q "Percolator version"
    "$DOTNET_EXE" --list-runtimes | grep -q "Microsoft.NETCore.App 8."
    test -f "$FLASHLFQ_EXE"
    test -f "$THERMO_EXE"
    mono --version | head -1
' || fail "tools missing in the worker image"

echo "OK: web UI smoke test passed"
