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

    return { store, submitting, error, success, submit };
}
