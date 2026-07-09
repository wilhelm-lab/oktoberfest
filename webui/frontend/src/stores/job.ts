import { defineStore } from "pinia";
import { ref, reactive } from "vue";
import * as api from "@/api/client";

export interface MetaCatalog {
    models: { intensity: string[]; irt: string[] };
    searchEngines: string[];
    spectraTypes: string[];
    libraryFormats: string[];
    enzymes: string[];
    tags: string[];
}

export interface UploadEntry {
    file: File;
    progress: number;
    status: "pending" | "uploading" | "done" | "error";
    error?: string;
}

export const useJobStore = defineStore("job", () => {
    const meta = reactive<MetaCatalog>({
        models: { intensity: [], irt: [] },
        searchEngines: [],
        spectraTypes: [],
        libraryFormats: [],
        enzymes: [],
        tags: [],
    });
    const metaLoaded = ref(false);

    // Per-role upload lists
    const uploads = reactive<Record<string, UploadEntry[]>>({});
    const isUploading = ref(false);

    // Current job being viewed
    const currentJob = ref<Record<string, unknown> | null>(null);

    // Job being drafted (created but not submitted)
    const activeDraftJobId = ref<string | null>(null);
    const activeDraftJobType = ref<string | null>(null);

    async function fetchMeta() {
        if (metaLoaded.value) return;
        const [models, engines, spectra, formats, enzymes, tags] =
            await Promise.all([
                api.getModels(),
                api.getSearchEngines(),
                api.getSpectraTypes(),
                api.getLibraryFormats(),
                api.getEnzymes(),
                api.getTags(),
            ]);
        meta.models = models;
        meta.searchEngines = engines;
        meta.spectraTypes = spectra;
        meta.libraryFormats = formats;
        meta.enzymes = enzymes;
        meta.tags = tags;
        metaLoaded.value = true;
    }

    async function fetchDefaults(jobType: string) {
        return api.getDefaults(jobType);
    }

    function clearUploads() {
        Object.keys(uploads).forEach((k) => delete uploads[k]);
        activeDraftJobId.value = null;
        activeDraftJobType.value = null;
    }

    function addUpload(role: string, file: File) {
        if (!uploads[role]) uploads[role] = [];
        uploads[role].push({ file, progress: 0, status: "pending" });
    }

    async function ensureJob(jobType: string): Promise<string> {
        if (activeDraftJobId.value && activeDraftJobType.value === jobType) {
            return activeDraftJobId.value;
        }
        // Create a job with empty config, since all defaults are valid
        const { job_id } = await api.createJob(jobType, {});
        activeDraftJobId.value = job_id;
        activeDraftJobType.value = jobType;
        return job_id;
    }

    async function uploadRole(jobType: string, role: string): Promise<void> {
        const jobId = await ensureJob(jobType);
        const entries = uploads[role];
        if (!entries || entries.length === 0) return;

        isUploading.value = true;
        try {
            for (const entry of entries) {
                if (entry.status === "done") continue;
                entry.status = "uploading";
                try {
                    await api.uploadFile(jobId, role, entry.file, (e) => {
                        if (e.total)
                            entry.progress = Math.round((e.loaded / e.total) * 100);
                    });
                    entry.status = "done";
                    entry.progress = 100;
                } catch (err: unknown) {
                    entry.status = "error";
                    entry.error = err instanceof Error ? err.message : String(err);
                    throw err;
                }
            }
        } finally {
            isUploading.value = false;
        }
    }

    async function submitDraftJob(config: Record<string, unknown>): Promise<string> {
        if (!activeDraftJobId.value) {
            throw new Error("No active job. Please upload files first.");
        }
        const jobId = activeDraftJobId.value;
        await api.submitJob(jobId, config);
        clearUploads();
        return jobId;
    }

    async function fetchStatus(jobId: string) {
        currentJob.value = await api.getJob(jobId);
        return currentJob.value;
    }

    return {
        meta,
        metaLoaded,
        uploads,
        isUploading,
        currentJob,
        activeDraftJobId,
        activeDraftJobType,
        fetchMeta,
        fetchDefaults,
        clearUploads,
        addUpload,
        ensureJob,
        uploadRole,
        submitDraftJob,
        fetchStatus,
    };
});
