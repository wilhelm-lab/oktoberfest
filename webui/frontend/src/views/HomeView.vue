<script setup lang="ts">
import { ref } from "vue";
import { useRouter } from "vue-router";

const router = useRouter();
const checkJobId = ref("");

function goToJob() {
    if (checkJobId.value.trim()) {
        router.push({
            name: "job",
            params: { jobId: checkJobId.value.trim() },
        });
    }
}

const jobTypes = [
    {
        name: "Rescoring",
        route: "rescoring",
        icon: "mdi-filter-variant",
        color: "#0E2D4E",
        description:
            "Rescore PSMs using deep learning-predicted spectra to improve FDR.",
    },
    {
        name: "CE Calibration",
        route: "ce",
        icon: "mdi-tune",
        color: "#1565C0",
        description:
            "Estimate the optimal normalized collision energy (NCE) for your dataset.",
    },
    {
        name: "Spectral Library",
        route: "speclib",
        icon: "mdi-library",
        color: "#2E7D32",
        description:
            "Build a predicted spectral library from a FASTA or peptide list.",
    },
];
</script>

<template>
    <v-container class="py-8" max-width="960">
        <div class="text-center mb-8">
            <h1 class="text-h4 font-weight-bold text-primary mb-2">
                Oktoberfest
            </h1>
            <p class="text-body-1 text-grey-darken-1">
                Open-source spectral library generation and rescoring pipeline.
                Submit and monitor jobs through your browser.
            </p>
        </div>

        <v-row class="mb-8">
            <v-col v-for="jt in jobTypes" :key="jt.name" cols="12" md="4">
                <v-card
                    :color="jt.color"
                    class="text-white pa-2"
                    hover
                    style="cursor: pointer"
                    @click="router.push({ name: jt.route })"
                >
                    <v-card-text class="text-center pa-6">
                        <v-icon :size="48" color="white" class="mb-3">{{
                            jt.icon
                        }}</v-icon>
                        <div class="text-h6 mb-2">{{ jt.name }}</div>
                        <div class="text-body-2 opacity-80">
                            {{ jt.description }}
                        </div>
                    </v-card-text>
                    <v-card-actions class="justify-center pb-4">
                        <v-btn variant="outlined" color="white">Start →</v-btn>
                    </v-card-actions>
                </v-card>
            </v-col>
        </v-row>

        <v-card variant="outlined" class="mb-4">
            <v-card-title class="text-body-1 font-weight-medium"
                >Check job status</v-card-title
            >
            <v-card-text>
                <div class="d-flex ga-3">
                    <v-text-field
                        v-model="checkJobId"
                        label="Job ID"
                        placeholder="Paste your Job ID here"
                        density="compact"
                        variant="outlined"
                        hide-details
                        @keyup.enter="goToJob"
                    />
                    <v-btn
                        color="primary"
                        :disabled="!checkJobId.trim()"
                        @click="goToJob"
                    >
                        Check status
                    </v-btn>
                </div>
            </v-card-text>
        </v-card>

        <div class="text-center">
            <v-btn variant="text" color="primary" :to="{ name: 'jobs' }">
                View recent jobs
            </v-btn>
        </div>
    </v-container>
</template>
