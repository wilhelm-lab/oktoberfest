<script setup lang="ts">
import { reactive, computed } from "vue";
import FileUpload from "@/components/FileUpload.vue";
import RangeField from "@/components/RangeField.vue";
import ModelsSection from "@/components/ModelsSection.vue";
import ConfigSummary from "@/components/ConfigSummary.vue";
import { useJobForm } from "@/composables/useJobForm";

const { store, submitting, error, submit } = useJobForm(
    "SpectralLibraryGeneration"
);

const form = reactive({
    library_input_type: "fasta" as "fasta" | "peptides" | "",
    intensity: "Prosit_2020_intensity_HCD",
    irt: "Prosit_2019_irt",
    prediction_server: "koina.wilhelmlab.org:443",
    ssl: true,
    numThreads: 1,
    tag: "",
    instrument_type: "QE",
    // spectralLibraryOptions
    fragmentation: "HCD",
    collisionEnergy: 30,
    precursorCharge: [2, 3] as [number, number],
    minIntensity: 0.0005,
    nrOx: 1,
    batchsize: 10000,
    format: "msp",
    // fastaDigestOptions
    digestion: "full",
    missedCleavages: 2,
    minLength: 7,
    maxLength: 60,
    enzyme: "trypsin",
    specialAas: "KR",
    db: "concat",
});

const fileAccept = computed(() =>
    form.library_input_type === "peptides"
        ? [".csv", ".tsv", ".txt"]
        : [".fasta", ".fa", ".faa"]
);

const uploadRole = computed(() =>
    form.library_input_type === "peptides" ? "peptides" : "fasta"
);

const libraryFiles = computed(
    () => store.uploads["fasta"] ?? store.uploads["peptides"] ?? []
);

const isValid = computed(
    () =>
        libraryFiles.value.some((f) => f.status !== "error") &&
        form.precursorCharge[0] < form.precursorCharge[1]
);

function buildConfig(): Record<string, unknown> {
    const cfg: Record<string, unknown> = {
        inputs: {
            library_input_type: form.library_input_type,
            instrument_type: form.instrument_type,
        },
        models: { intensity: form.intensity, irt: form.irt },
        prediction_server: form.prediction_server,
        ssl: form.ssl,
        numThreads: form.numThreads,
        tag: form.tag,
        spectralLibraryOptions: {
            fragmentation: form.fragmentation,
            collisionEnergy: form.collisionEnergy,
            precursorCharge: form.precursorCharge,
            minIntensity: form.minIntensity,
            nrOx: form.nrOx,
            batchsize: form.batchsize,
            format: form.format,
        },
    };
    if (form.library_input_type === "fasta") {
        cfg["fastaDigestOptions"] = {
            fragmentation: form.fragmentation,
            digestion: form.digestion,
            missedCleavages: form.missedCleavages,
            minLength: form.minLength,
            maxLength: form.maxLength,
            enzyme: form.enzyme,
            specialAas: form.specialAas,
            db: form.db,
        };
    }
    return cfg;
}

function onFilesAdded(role: string, files: File[]) {
    files.forEach((f) => store.addUpload(role, f));
}

async function handleSubmit() {
    await submit(buildConfig());
}
</script>

<template>
    <v-container max-width="900" class="py-6">
        <h2 class="text-h5 font-weight-bold text-primary mb-6">
            Spectral Library Generation
        </h2>

        <v-alert
            v-if="error"
            type="error"
            class="mb-4"
            closable
            @click:close="error = ''"
            >{{ error }}</v-alert
        >

        <!-- Inputs -->
        <v-card variant="outlined" class="mb-4">
            <v-card-title class="text-body-1 font-weight-medium"
                >Inputs & Files</v-card-title
            >
            <v-card-text>
                <v-row>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.library_input_type"
                            :items="[
                                {
                                    title: 'FASTA (in-silico digest)',
                                    value: 'fasta',
                                },
                                { title: 'Peptides CSV', value: 'peptides' },
                                { title: 'Internal format', value: '' },
                            ]"
                            label="Library source *"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12">
                        <div class="text-body-2 mb-1 font-weight-medium">
                            {{
                                form.library_input_type === "peptides"
                                    ? "Peptides CSV"
                                    : "FASTA file"
                            }}
                            *
                        </div>
                        <FileUpload
                            :role="uploadRole"
                            :accept="fileAccept"
                            :multiple="false"
                            @add="onFilesAdded(uploadRole, $event)"
                        />
                    </v-col>
                </v-row>
            </v-card-text>
        </v-card>

        <!-- Models -->
        <ModelsSection
            v-model:intensity-model="form.intensity"
            v-model:irt-model="form.irt"
            v-model:prediction-server="form.prediction_server"
            v-model:ssl="form.ssl"
        />

        <!-- Library options -->
        <v-card variant="outlined" class="mb-4">
            <v-card-title class="text-body-1 font-weight-medium"
                >Library Options</v-card-title
            >
            <v-card-text>
                <v-row>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.fragmentation"
                            :items="['HCD', 'CID']"
                            label="Fragmentation"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.collisionEnergy"
                            label="Collision energy"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.format"
                            :items="store.meta.libraryFormats"
                            label="Output format"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12">
                        <RangeField
                            v-model="form.precursorCharge"
                            label="Precursor charge range (min, max)"
                            :min="1"
                            :max="10"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.minIntensity"
                            label="Min relative intensity"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.nrOx"
                            label="Max oxidations (M)"
                            type="number"
                            density="compact"
                            variant="outlined"
                            :min="0"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.batchsize"
                            label="Batch size"
                            type="number"
                            density="compact"
                            variant="outlined"
                            :min="1"
                        />
                    </v-col>
                </v-row>
            </v-card-text>
        </v-card>

        <!-- Digestion (shown only for FASTA) -->
        <v-card
            v-if="form.library_input_type === 'fasta'"
            variant="outlined"
            class="mb-4"
        >
            <v-card-title class="text-body-1 font-weight-medium"
                >In-silico Digestion</v-card-title
            >
            <v-card-text>
                <v-row>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.enzyme"
                            :items="store.meta.enzymes"
                            label="Enzyme"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.digestion"
                            :items="['full', 'semi', 'none']"
                            label="Digestion mode"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.missedCleavages"
                            label="Missed cleavages"
                            type="number"
                            density="compact"
                            variant="outlined"
                            :min="0"
                        />
                    </v-col>
                    <v-col cols="12" md="3">
                        <v-text-field
                            v-model.number="form.minLength"
                            label="Min peptide length"
                            type="number"
                            density="compact"
                            variant="outlined"
                            :min="1"
                        />
                    </v-col>
                    <v-col cols="12" md="3">
                        <v-text-field
                            v-model.number="form.maxLength"
                            label="Max peptide length"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="3">
                        <v-text-field
                            v-model="form.specialAas"
                            label="Special AAs (decoys)"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="3">
                        <v-select
                            v-model="form.db"
                            :items="['target', 'decoy', 'concat']"
                            label="Database"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                </v-row>
            </v-card-text>
        </v-card>

        <!-- Advanced -->
        <v-expansion-panels class="mb-4">
            <v-expansion-panel>
                <v-expansion-panel-title
                    >Advanced settings</v-expansion-panel-title
                >
                <v-expansion-panel-text>
                    <v-row>
                        <v-col cols="12" md="4"
                            ><v-text-field
                                v-model.number="form.numThreads"
                                label="Parallel threads"
                                type="number"
                                density="compact"
                                variant="outlined"
                                :min="1"
                        /></v-col>
                        <v-col cols="12" md="4"
                            ><v-select
                                v-model="form.tag"
                                :items="store.meta.tags"
                                label="Isobaric tag"
                                density="compact"
                                variant="outlined"
                        /></v-col>
                        <v-col cols="12" md="4"
                            ><v-select
                                v-model="form.instrument_type"
                                :items="['QE', 'LUMOS', 'TIMSTOF', 'SCIEXTOF']"
                                label="Instrument type"
                                density="compact"
                                variant="outlined"
                        /></v-col>
                    </v-row>
                </v-expansion-panel-text>
            </v-expansion-panel>
        </v-expansion-panels>

        <ConfigSummary :config="buildConfig()" class="mb-4" />

        <div class="d-flex justify-end">
            <v-btn
                color="primary"
                size="large"
                :loading="submitting"
                :disabled="!isValid || submitting"
                @click="handleSubmit"
            >
                Submit job
            </v-btn>
        </div>
    </v-container>
</template>
