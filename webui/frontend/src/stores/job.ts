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

    // Current job being viewed
    const currentJob = ref<Record<string, unknown> | null>(null);

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
    }

    function addUpload(role: string, file: File) {
        if (!uploads[role]) uploads[role] = [];
        uploads[role].push({ file, progress: 0, status: "pending" });
    }

    async function uploadFiles(jobId: string): Promise<void> {
        for (const [role, entries] of Object.entries(uploads)) {
            for (const entry of entries) {
                entry.status = "uploading";
                try {
                    await api.uploadFile(jobId, role, entry.file, (e) => {
                        if (e.total)
                            entry.progress = Math.round(
                                (e.loaded / e.total) * 100
                            );
                    });
                    entry.status = "done";
                    entry.progress = 100;
                } catch (err: unknown) {
                    entry.status = "error";
                    entry.error =
                        err instanceof Error ? err.message : String(err);
                    throw err;
                }
            }
        }
    }

    async function createAndSubmit(
        jobType: string,
        config: Record<string, unknown>
    ): Promise<string> {
        const { job_id } = await api.createJob(jobType, config);
        await uploadFiles(job_id);
        await api.submitJob(job_id);
        clearUploads();
        return job_id;
    }

    async function fetchStatus(jobId: string) {
        currentJob.value = await api.getJob(jobId);
        return currentJob.value;
    }

    return {
        meta,
        metaLoaded,
        uploads,
        currentJob,
        fetchMeta,
        fetchDefaults,
        clearUploads,
        addUpload,
        uploadFiles,
        createAndSubmit,
        fetchStatus,
    };
});
