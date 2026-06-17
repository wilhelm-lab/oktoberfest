<script setup lang="ts">
import { useJobStore } from "@/stores/job";

const store = useJobStore();

defineProps<{
    intensityModel: string;
    irtModel: string;
    predictionServer: string;
    ssl: boolean;
}>();

const emit = defineEmits<{
    (e: "update:intensityModel", v: string): void;
    (e: "update:irtModel", v: string): void;
    (e: "update:predictionServer", v: string): void;
    (e: "update:ssl", v: boolean): void;
}>();
</script>

<template>
    <v-card variant="outlined" class="mb-4">
        <v-card-title class="text-body-1 font-weight-medium"
            >Models & Prediction</v-card-title
        >
        <v-card-text>
            <v-row>
                <v-col cols="12" md="6">
                    <v-combobox
                        :model-value="intensityModel"
                        :items="store.meta.models.intensity"
                        label="Intensity model"
                        density="compact"
                        variant="outlined"
                        @update:model-value="
                            emit('update:intensityModel', $event)
                        "
                    />
                </v-col>
                <v-col cols="12" md="6">
                    <v-combobox
                        :model-value="irtModel"
                        :items="['', ...store.meta.models.irt]"
                        label="iRT model (leave blank to disable)"
                        density="compact"
                        variant="outlined"
                        @update:model-value="
                            emit('update:irtModel', $event ?? '')
                        "
                    />
                </v-col>
                <v-col cols="12" md="8">
                    <v-text-field
                        :model-value="predictionServer"
                        label="Prediction server (host:port)"
                        density="compact"
                        variant="outlined"
                        @update:model-value="
                            emit('update:predictionServer', $event)
                        "
                    />
                </v-col>
                <v-col cols="12" md="4" class="d-flex align-center">
                    <v-switch
                        :model-value="ssl"
                        label="Use SSL"
                        color="primary"
                        hide-details
                        @update:model-value="emit('update:ssl', !!$event)"
                    />
                </v-col>
            </v-row>
        </v-card-text>
    </v-card>
</template>
