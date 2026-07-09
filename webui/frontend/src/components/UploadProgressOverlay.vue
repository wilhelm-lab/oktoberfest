<script setup lang="ts">
import { computed } from "vue";
import { useJobStore } from "@/stores/job";

const store = useJobStore();

const allUploads = computed(() => {
    return Object.values(store.uploads).flat();
});

const activeUpload = computed(() => {
    return allUploads.value.find((u) => u.status === "uploading");
});

const uploadProgress = computed(() => {
    return activeUpload.value ? activeUpload.value.progress : 0;
});
</script>

<template>
    <v-dialog :model-value="store.isUploading" persistent max-width="500">
        <v-card class="pa-4">
            <v-card-title class="text-h6 pb-2 text-primary font-weight-bold">
                Uploading files...
            </v-card-title>
            <v-card-text>
                <div class="mb-2 text-body-2">
                    <span v-if="activeUpload">
                        Currently uploading: <strong>{{ activeUpload.file.name }}</strong>
                    </span>
                    <span v-else>Processing...</span>
                </div>
                <v-progress-linear
                    :model-value="uploadProgress"
                    color="primary"
                    height="20"
                    striped
                    animated
                    rounded
                >
                    <template #default="{ value }">
                        <strong>{{ Math.ceil(value) }}%</strong>
                    </template>
                </v-progress-linear>
                
                <div class="mt-4 text-caption text-grey">
                    Please do not close this window until the upload is complete.
                </div>
            </v-card-text>
        </v-card>
    </v-dialog>
</template>
