<script setup lang="ts">
import { ref, computed, onMounted, onUnmounted } from "vue";
import { useRoute } from "vue-router";
import { getJob, getJobLog, resultsUrl } from "@/api/client";

const route = useRoute();
const jobId = computed(() => route.params.jobId as string);

const job = ref<Record<string, unknown> | null>(null);
const log = ref("");
const loading = ref(true);
const loadingLog = ref(false);
const error = ref("");
const showLog = ref(false);
const copied = ref(false);
const checkOtherId = ref("");

const TERMINAL = new Set(["SUCCEEDED", "FAILED", "CANCELLED"]);

let pollInterval: ReturnType<typeof setInterval> | null = null;

async function fetchJob() {
    try {
        job.value = await getJob(jobId.value);
        if (job.value && TERMINAL.has(job.value.status as string)) {
            stopPolling();
            if (job.value.status === "FAILED") {
                await fetchLog();
            }
        }
    } catch (e: unknown) {
        error.value = e instanceof Error ? e.message : String(e);
        stopPolling();
    } finally {
        loading.value = false;
    }
}

async function fetchLog() {
    loadingLog.value = true;
    try {
        const data = await getJobLog(jobId.value);
        log.value = data.log ?? "";
    } finally {
        loadingLog.value = false;
    }
}

function stopPolling() {
    if (pollInterval) {
        clearInterval(pollInterval);
        pollInterval = null;
    }
}

function startPolling() {
    pollInterval = setInterval(fetchJob, 3500);
}

onMounted(() => {
    fetchJob();
    startPolling();
});
onUnmounted(stopPolling);

function copyJobId() {
    navigator.clipboard.writeText(jobId.value);
    copied.value = true;
    setTimeout(() => (copied.value = false), 2000);
}

function download() {
    window.open(resultsUrl(jobId.value), "_blank");
}

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
    <v-container max-width="800" class="py-6">
        <v-btn
            variant="text"
            :to="{ name: 'home' }"
            prepend-icon="mdi-arrow-left"
            class="mb-4"
        >
            Home
        </v-btn>

        <v-skeleton-loader v-if="loading" type="card" />

        <v-alert v-else-if="error" type="error" class="mb-4">{{
            error
        }}</v-alert>

        <template v-else-if="job">
            <h2 class="text-h5 font-weight-bold text-primary mb-2">
                {{ job.job_type }} Job
            </h2>

            <!-- Job ID + copy -->
            <div class="d-flex align-center ga-2 mb-4">
                <v-chip label color="primary" size="small" class="font-mono">{{
                    jobId
                }}</v-chip>
                <v-btn
                    icon="mdi-content-copy"
                    size="x-small"
                    variant="text"
                    @click="copyJobId"
                    :color="copied ? 'success' : undefined"
                />
                <span v-if="copied" class="text-caption text-success"
                    >Copied!</span
                >
            </div>

            <!-- Status chip -->
            <v-chip
                :color="statusColor[job.status as string] ?? 'grey'"
                class="mb-4"
                label
                size="large"
            >
                <v-icon start>
                    {{
                        job.status === "RUNNING"
                            ? "mdi-loading mdi-spin"
                            : job.status === "SUCCEEDED"
                            ? "mdi-check-circle"
                            : job.status === "FAILED"
                            ? "mdi-alert-circle"
                            : "mdi-clock-outline"
                    }}
                </v-icon>
                {{ job.status }}
            </v-chip>

            <!-- Timestamps -->
            <v-table density="compact" class="mb-4">
                <tbody>
                    <tr>
                        <td class="text-grey">Created</td>
                        <td>{{ fmtDate(job.created_at as string) }}</td>
                    </tr>
                    <tr>
                        <td class="text-grey">Started</td>
                        <td>{{ fmtDate(job.started_at as string) }}</td>
                    </tr>
                    <tr>
                        <td class="text-grey">Finished</td>
                        <td>{{ fmtDate(job.finished_at as string) }}</td>
                    </tr>
                </tbody>
            </v-table>

            <!-- QUEUED / RUNNING message -->
            <v-alert
                v-if="job.status === 'QUEUED' || job.status === 'RUNNING'"
                type="info"
                class="mb-4"
            >
                <v-progress-circular
                    v-if="job.status === 'RUNNING'"
                    indeterminate
                    size="18"
                    class="mr-2"
                />
                In progress. Large jobs can take hours. Save this URL or your
                Job ID to check back later.
            </v-alert>

            <!-- SUCCEEDED: download button -->
            <v-btn
                v-if="job.status === 'SUCCEEDED'"
                color="success"
                size="large"
                prepend-icon="mdi-download"
                class="mb-4"
                @click="download"
            >
                Download results.zip
            </v-btn>

            <!-- FAILED: error + log -->
            <template v-if="job.status === 'FAILED'">
                <v-alert type="error" class="mb-4">
                    <strong>Job failed:</strong> {{ job.error }}
                </v-alert>
                <v-expansion-panels v-model="showLog" class="mb-4">
                    <v-expansion-panel>
                        <v-expansion-panel-title>
                            <v-icon start>mdi-text-box</v-icon> View captured
                            log
                        </v-expansion-panel-title>
                        <v-expansion-panel-text>
                            <v-progress-circular
                                v-if="loadingLog"
                                indeterminate
                                size="24"
                            />
                            <pre
                                v-else
                                style="
                                    white-space: pre-wrap;
                                    font-size: 12px;
                                    max-height: 400px;
                                    overflow: auto;
                                "
                                >{{ log || "(no log available)" }}</pre
                            >
                        </v-expansion-panel-text>
                    </v-expansion-panel>
                </v-expansion-panels>
            </template>
        </template>

        <!-- Check another job -->
        <v-divider class="my-6" />
        <v-card variant="outlined">
            <v-card-title class="text-body-2 font-weight-medium"
                >Check another job</v-card-title
            >
            <v-card-text>
                <div class="d-flex ga-2">
                    <v-text-field
                        v-model="checkOtherId"
                        label="Job ID"
                        density="compact"
                        variant="outlined"
                        hide-details
                        @keyup.enter="
                            checkOtherId &&
                                $router.push({
                                    name: 'job',
                                    params: { jobId: checkOtherId },
                                })
                        "
                    />
                    <v-btn
                        color="primary"
                        :disabled="!checkOtherId.trim()"
                        @click="
                            $router.push({
                                name: 'job',
                                params: { jobId: checkOtherId },
                            })
                        "
                        >Go</v-btn
                    >
                </div>
            </v-card-text>
        </v-card>
    </v-container>
</template>
