<script setup lang="ts">
import { ref, onMounted } from "vue";
import { useRouter } from "vue-router";
import { getJobs } from "@/api/client";

const router = useRouter();
const jobs = ref<Record<string, unknown>[]>([]);
const loading = ref(true);
const error = ref("");

onMounted(async () => {
    try {
        jobs.value = await getJobs(50, 0);
    } catch (e: unknown) {
        error.value = e instanceof Error ? e.message : String(e);
    } finally {
        loading.value = false;
    }
});

const statusColor: Record<string, string> = {
    CREATED: "grey",
    QUEUED: "blue",
    RUNNING: "orange",
    SUCCEEDED: "green",
    FAILED: "red",
    CANCELLED: "grey",
};

function fmtDate(d: string | null | undefined) {
    if (!d) return "—";
    return new Date(d + "Z").toLocaleString();
}
</script>

<template>
    <v-container max-width="960" class="py-6">
        <h2 class="text-h5 font-weight-bold text-primary mb-4">Recent Jobs</h2>

        <v-alert v-if="error" type="error" class="mb-4">{{ error }}</v-alert>
        <v-skeleton-loader v-if="loading" type="table" />

        <v-table v-else hover>
            <thead>
                <tr>
                    <th>Job ID</th>
                    <th>Type</th>
                    <th>Status</th>
                    <th>Created</th>
                    <th>Finished</th>
                </tr>
            </thead>
            <tbody>
                <tr
                    v-for="job in jobs"
                    :key="job.job_id as string"
                    style="cursor: pointer"
                    @click="
                        router.push({
                            name: 'job',
                            params: { jobId: String(job.job_id) },
                        })
                    "
                >
                    <td class="font-mono text-caption">
                        {{ (job.job_id as string).slice(0, 8) }}…
                    </td>
                    <td>{{ job.job_type }}</td>
                    <td>
                        <v-chip
                            :color="statusColor[job.status as string] ?? 'grey'"
                            size="x-small"
                            label
                        >
                            {{ job.status }}
                        </v-chip>
                    </td>
                    <td class="text-caption">
                        {{ fmtDate(job.created_at as string) }}
                    </td>
                    <td class="text-caption">
                        {{ fmtDate(job.finished_at as string) }}
                    </td>
                </tr>
            </tbody>
        </v-table>

        <div v-if="!loading && !jobs.length" class="text-center text-grey py-8">
            No jobs yet.
            <router-link to="/">Submit your first job</router-link>.
        </div>
    </v-container>
</template>
