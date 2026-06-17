import axios, { type AxiosProgressEvent } from "axios";

const api = axios.create({
    baseURL: import.meta.env.VITE_API_BASE ?? "/api/v1",
    timeout: 30000,
});

// Centralized error handling
api.interceptors.response.use(
    (r) => r,
    (err) => {
        const msg =
            err.response?.data?.detail ?? err.message ?? "Network error";
        console.error("[API]", msg);
        return Promise.reject(
            new Error(typeof msg === "string" ? msg : JSON.stringify(msg))
        );
    }
);

// ── Job endpoints ──────────────────────────────────────────────────────────────

export function createJob(jobType: string, config: Record<string, unknown>) {
    return api.post("/jobs", { job_type: jobType, config }).then((r) => r.data);
}

export function uploadFile(
    jobId: string,
    role: string,
    file: File,
    onProgress?: (e: AxiosProgressEvent) => void
) {
    const fd = new FormData();
    fd.append("role", role);
    fd.append("file", file);
    return api
        .post(`/jobs/${jobId}/files`, fd, {
            headers: { "Content-Type": "multipart/form-data" },
            onUploadProgress: onProgress,
        })
        .then((r) => r.data);
}

export function submitJob(jobId: string) {
    return api.post(`/jobs/${jobId}/submit`).then((r) => r.data);
}

export function getJob(jobId: string) {
    return api.get(`/jobs/${jobId}`).then((r) => r.data);
}

export function getJobLog(jobId: string) {
    return api.get(`/jobs/${jobId}/log`).then((r) => r.data);
}

export function getJobs(limit = 50, offset = 0) {
    return api.get("/jobs", { params: { limit, offset } }).then((r) => r.data);
}

export function cancelJob(jobId: string) {
    return api.post(`/jobs/${jobId}/cancel`).then((r) => r.data);
}

export function resultsUrl(jobId: string): string {
    const base = import.meta.env.VITE_API_BASE ?? "/api/v1";
    return `${base}/jobs/${jobId}/results`;
}

// ── Meta endpoints ─────────────────────────────────────────────────────────────

export function getModels() {
    return api.get("/meta/models").then((r) => r.data);
}

export function getSearchEngines() {
    return api.get("/meta/search-engines").then((r) => r.data);
}

export function getSpectraTypes() {
    return api.get("/meta/spectra-types").then((r) => r.data);
}

export function getLibraryFormats() {
    return api.get("/meta/library-formats").then((r) => r.data);
}

export function getEnzymes() {
    return api.get("/meta/enzymes").then((r) => r.data);
}

export function getTags() {
    return api.get("/meta/tags").then((r) => r.data);
}

export function getDefaults(jobType: string) {
    return api.get(`/meta/defaults/${jobType}`).then((r) => r.data);
}

export function getHealth() {
    return api.get("/health").then((r) => r.data);
}
