import { ref, reactive, onMounted } from "vue";
import { useRouter } from "vue-router";
import { useJobStore } from "@/stores/job";

export function useJobForm(jobType: string) {
    const store = useJobStore();
    const router = useRouter();
    const submitting = ref(false);
    const error = ref("");
    const success = ref("");

    onMounted(() => store.fetchMeta());

    async function submit(config: Record<string, unknown>) {
        error.value = "";
        submitting.value = true;
        try {
            const jobId = await store.submitDraftJob(config);
            router.push({ name: "job", params: { jobId } });
        } catch (err: unknown) {
            error.value = err instanceof Error ? err.message : String(err);
        } finally {
            submitting.value = false;
        }
    }

    // Pre-fill fields whose default is configured on the server (e.g. PREDICTION_SERVER)
    async function loadServerDefaults(form: { prediction_server: string }) {
        try {
            const defaults = await store.fetchDefaults(jobType);
            if (defaults?.prediction_server) form.prediction_server = defaults.prediction_server;
        } catch {
            // keep the form's built-in default
        }
    }

    return { store, submitting, error, success, submit, loadServerDefaults };
}
