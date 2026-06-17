<script setup lang="ts">
import { reactive, computed } from "vue";
import FileUpload from "@/components/FileUpload.vue";
import RangeField from "@/components/RangeField.vue";
import ModelsSection from "@/components/ModelsSection.vue";
import ConfigSummary from "@/components/ConfigSummary.vue";
import { useJobForm } from "@/composables/useJobForm";

const { store, submitting, error, submit } = useJobForm("Rescoring");

const form = reactive({
    // Inputs
    search_results_type: "Maxquant",
    spectra_type: "mzml",
    // Models
    intensity: "Prosit_2020_intensity_HCD",
    irt: "Prosit_2019_irt",
    prediction_server: "koina.wilhelmlab.org:443",
    ssl: true,
    // Job options
    fdr_estimation_method: "mokapot",
    regressionMethod: "spline",
    add_feature_cols: "none",
    massTolerance: 20,
    unitMassTolerance: "ppm",
    ce_range: [19, 50] as [number, number],
    use_ransac_model: false,
    thermoExe: "ThermoRawFileParser.exe",
    // Quantification
    quantification: false,
    // Advanced
    numThreads: 1,
    tag: "",
    instrument_type: "QE",
    ion_types_list: ["b", "y"],
});

const searchResultFiles = computed(() => store.uploads["search_results"] ?? []);
const spectraFiles = computed(() => store.uploads["spectra"] ?? []);
const fastaFiles = computed(() => store.uploads["fasta"] ?? []);

const showThermoExe = computed(() => form.spectra_type === "raw");

const isValid = computed(() => {
    return (
        searchResultFiles.value.some((f) => f.status !== "error") &&
        spectraFiles.value.some((f) => f.status !== "error") &&
        form.ce_range[0] < form.ce_range[1]
    );
});

function buildConfig(): Record<string, unknown> {
    const cfg: Record<string, unknown> = {
        inputs: {
            search_results_type: form.search_results_type,
            spectra_type: form.spectra_type,
            instrument_type: form.instrument_type,
        },
        models: { intensity: form.intensity, irt: form.irt },
        prediction_server: form.prediction_server,
        ssl: form.ssl,
        numThreads: form.numThreads,
        tag: form.tag,
        fdr_estimation_method: form.fdr_estimation_method,
        regressionMethod: form.regressionMethod,
        add_feature_cols: form.add_feature_cols,
        massTolerance: form.massTolerance,
        unitMassTolerance: form.unitMassTolerance,
        ion_types: [...form.ion_types_list].sort().reverse().join(""),
        ce_alignment_options: {
            ce_range: form.ce_range,
            use_ransac_model: form.use_ransac_model,
        },
        quantification: form.quantification,
    };
    if (showThermoExe.value) cfg["thermoExe"] = form.thermoExe;
    if (form.quantification) {
        cfg["fastaDigestOptions"] = {
            digestion: "full",
            missedCleavages: 2,
            minLength: 7,
            maxLength: 60,
            enzyme: "trypsin",
            specialAas: "KR",
            db: "concat",
            fragmentation: "HCD",
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
        <h2 class="text-h5 font-weight-bold text-primary mb-6">Rescoring</h2>

        <v-alert
            v-if="error"
            type="error"
            class="mb-4"
            closable
            @click:close="error = ''"
        >
            {{ error }}
        </v-alert>

        <!-- Inputs & Files -->
        <v-card variant="outlined" class="mb-4">
            <v-card-title class="text-body-1 font-weight-medium"
                >Inputs & Files</v-card-title
            >
            <v-card-text>
                <v-row>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.search_results_type"
                            :items="store.meta.searchEngines"
                            label="Search engine"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.spectra_type"
                            :items="store.meta.spectraTypes"
                            label="Spectra type"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12">
                        <div class="text-body-2 mb-1 font-weight-medium">
                            Search results file(s) *
                        </div>
                        <FileUpload
                            role="search_results"
                            :accept="[
                                '.txt',
                                '.tsv',
                                '.csv',
                                '.pepxml',
                                '.mzid',
                                '.xml',
                            ]"
                            :multiple="true"
                            @add="onFilesAdded('search_results', $event)"
                        />
                    </v-col>
                    <v-col cols="12">
                        <div class="text-body-2 mb-1 font-weight-medium">
                            Spectra file(s) *
                        </div>
                        <FileUpload
                            role="spectra"
                            :accept="[
                                '.raw',
                                '.mzML',
                                '.mzml',
                                '.hdf',
                                '.hdf5',
                                '.d',
                            ]"
                            :multiple="true"
                            @add="onFilesAdded('spectra', $event)"
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

        <!-- Job options -->
        <v-card variant="outlined" class="mb-4">
            <v-card-title class="text-body-1 font-weight-medium"
                >Job Options</v-card-title
            >
            <v-card-text>
                <v-row>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.fdr_estimation_method"
                            :items="['mokapot', 'percolator']"
                            label="FDR method"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.regressionMethod"
                            :items="['spline', 'lowess', 'logistic']"
                            label="iRT regression"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.add_feature_cols"
                            :items="['none', 'all']"
                            label="Extra feature columns"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.massTolerance"
                            label="Mass tolerance"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.unitMassTolerance"
                            :items="['ppm', 'da']"
                            label="Tolerance unit"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4" class="d-flex align-center">
                        <v-switch
                            v-model="form.use_ransac_model"
                            label="RANSAC (timsTOF)"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.ion_types_list"
                            :items="['a', 'b', 'c', 'x', 'y', 'z']"
                            label="Ion types for annotation"
                            density="compact"
                            variant="outlined"
                            multiple
                            chips
                        />
                    </v-col>
                    <v-col cols="12">
                        <RangeField
                            v-model="form.ce_range"
                            label="CE range (min, max)"
                            :min="1"
                            :max="100"
                        />
                    </v-col>
                    <v-col v-if="showThermoExe" cols="12">
                        <v-text-field
                            v-model="form.thermoExe"
                            label="ThermoRawFileParser path"
                            density="compact"
                            variant="outlined"
                            hint="Required for .raw spectra files"
                            persistent-hint
                        />
                    </v-col>
                </v-row>
            </v-card-text>
        </v-card>

        <!-- Quantification -->
        <v-card variant="outlined" class="mb-4">
            <v-card-title
                class="d-flex align-center text-body-1 font-weight-medium"
            >
                <span>Quantification</span>
                <v-switch
                    v-model="form.quantification"
                    color="primary"
                    hide-details
                    class="ml-4"
                />
            </v-card-title>
            <v-card-text v-if="form.quantification">
                <div class="text-body-2 mb-2 font-weight-medium">
                    FASTA file (for quant) *
                </div>
                <FileUpload
                    role="fasta"
                    :accept="['.fasta', '.fa', '.faa']"
                    :multiple="false"
                    @add="onFilesAdded('fasta', $event)"
                />
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
                        <v-col cols="12" md="4">
                            <v-text-field
                                v-model.number="form.numThreads"
                                label="Parallel threads"
                                type="number"
                                density="compact"
                                variant="outlined"
                                :min="1"
                            />
                        </v-col>
                        <v-col cols="12" md="4">
                            <v-select
                                v-model="form.tag"
                                :items="store.meta.tags"
                                label="Isobaric tag"
                                density="compact"
                                variant="outlined"
                            />
                        </v-col>
                        <v-col cols="12" md="4">
                            <v-select
                                v-model="form.instrument_type"
                                :items="['QE', 'LUMOS', 'TIMSTOF', 'SCIEXTOF']"
                                label="Instrument type"
                                density="compact"
                                variant="outlined"
                            />
                        </v-col>
                    </v-row>
                </v-expansion-panel-text>
            </v-expansion-panel>
        </v-expansion-panels>

        <!-- Config preview -->
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
