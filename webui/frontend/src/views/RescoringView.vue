<script setup lang="ts">
import { reactive, computed, onMounted, watch } from "vue";
import FileUpload from "@/components/FileUpload.vue";
import RangeField from "@/components/RangeField.vue";
import ModelsSection from "@/components/ModelsSection.vue";
import ConfigSummary from "@/components/ConfigSummary.vue";
import { useJobForm } from "@/composables/useJobForm";

const { store, submitting, error, submit, loadServerDefaults } = useJobForm("Rescoring");

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
    fdr_estimation_method: "percolator",
    regressionMethod: "spline",
    add_feature_cols: "none",
    massTolerance: 20,
    massToleranceAuto: false,
    unitMassTolerance: "ppm",
    ce_range: [19, 50] as [number, number],
    use_ransac_model: false,
    percolatorSubsetEnabled: true,
    percolatorSubsetMaxTrain: 500000,
    chimericPeakSubtraction: false,
    ptm_fdr_control: "auto" as "auto" | "on" | "off",
    // Quantification
    quantification: false,
    quant_precursor: "search_engine",
    quant_protein: "directlfq",
    quant_normalization: "stable_drift",
    quant_psmFdr: 0.01,
    quant_peptideFdr: 0.01,
    quant_proteinFdr: 0.01,
    quant_ppmTolerance: 10,
    quant_ppmFilterMode: "auto" as "auto" | "value" | "none",
    quant_ppmFilter: 5,
    quant_mbr: false,
    quant_mbrFdr: 0.01,
    digest_enzyme: "trypsin",
    digest_digestion: "full",
    digest_missedCleavages: 2,
    digest_minLength: 7,
    digest_maxLength: 60,
    digest_specialAas: "KR",
    // PTM localization
    ptm_localization: false,
    ...publishedLocalizationDefaults(),
    // Other Options
    tag: "",
    instrument_type: "QE",
    ion_types_list: ["b", "y"],
});
onMounted(() => loadServerDefaults(form));

// localization settings of the published benchmark runs (phospho_repro final model, docs/config.rst defaults)
function publishedLocalizationDefaults() {
    return {
        loc_ptms: ["21"] as string[],
        loc_residues: { "21": ["S", "T", "Y"] } as Record<string, string[]>,
        loc_neutral_loss: true,
        loc_decoy_residues: "A, G, L, V",
        loc_decoy_neutral_loss: false,
        loc_decoy_proxy: "",
        loc_decoy_selection: "residue",
        loc_flr_neighbor_split: false,
        loc_probability: "pep",
        loc_decoy_weighting: "count",
        loc_flr_groups: "none",
        loc_flr_estimator: "hybrid",
        loc_sdi_features: true,
        loc_ambiguity_check: true,
    };
}

function resetLocalizationDefaults() {
    Object.assign(form, publishedLocalizationDefaults());
}

const searchResultFiles = computed(() => store.uploads["search_results"] ?? []);
const spectraFiles = computed(() => store.uploads["spectra"] ?? []);
const fastaFiles = computed(() => store.uploads["fasta"] ?? []);

const AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY".split("");
const SEARCH_ENGINE_QUANT_ENGINES = ["Maxquant", "Sage"];

const ptmItems = computed(() =>
    Object.entries(store.meta.localizationPtms).map(([id, ptm]) => ({
        title: `${ptm.name} (UNIMOD:${id})`,
        value: id,
    }))
);

const precursorItems = computed(() =>
    (store.meta.quantification.precursor.length ? store.meta.quantification.precursor : ["search_engine", "flashlfq"]).map(
        (v) => ({
            title: v === "flashlfq" ? "FlashLFQ (MS1, any search engine)" : "Search engine (MaxQuant / Sage)",
            value: v,
            props: {
                disabled:
                    v === "flashlfq" &&
                    (!store.meta.quantification.flashlfq_available || form.ptm_localization),
            },
        })
    )
);

const decoyResidueList = computed(() =>
    form.loc_decoy_residues
        .toUpperCase()
        .split(/[\s,;]+/)
        .filter((aa) => aa.length > 0)
);

const possibleSites = computed(() => [
    ...new Set(form.loc_ptms.flatMap((id) => form.loc_residues[id] ?? [])),
]);

const isFlashlfq = computed(() => form.quant_precursor === "flashlfq");

// warnings shown above the submit button; any of them blocks submission
const quantProblem = computed(() => {
    if (!form.quantification) return "";
    if (isFlashlfq.value) {
        if (!store.meta.quantification.flashlfq_available)
            return "FlashLFQ is not configured on this server.";
        if (form.ptm_localization)
            return "FlashLFQ quantification cannot be combined with PTM localization.";
        if (!["raw", "mzml"].includes(form.spectra_type))
            return "FlashLFQ quantification needs raw or mzML spectra (MS1 scans).";
    } else if (!SEARCH_ENGINE_QUANT_ENGINES.includes(form.search_results_type)) {
        return "Search engine quantification supports MaxQuant (msms.txt + evidence.txt) and Sage (results.sage.tsv + lfq.tsv); use FlashLFQ for other search engines.";
    }
    return "";
});

const localizationProblem = computed(() => {
    if (!form.ptm_localization) return "";
    if (form.loc_ptms.length === 0) return "Select at least one PTM to localize.";
    for (const id of form.loc_ptms) {
        if (!(form.loc_residues[id] ?? []).length) return `Select the possible residues of UNIMOD:${id}.`;
    }
    const bad = decoyResidueList.value.filter((aa) => !AMINO_ACIDS.includes(aa));
    if (bad.length) return `Invalid decoy residues: ${bad.join(", ")}.`;
    const overlap = decoyResidueList.value.filter((aa) => possibleSites.value.includes(aa));
    if (overlap.length) return `Decoy residues ${overlap.join(", ")} can also carry a localized PTM.`;
    if (!form.intensity.includes("PTM")) return "PTM localization needs a PTM intensity model (e.g. Prosit_2024_intensity_PTMs_gl).";
    return "";
});

// localization needs a PTM model: switch to the PTM models when it is turned on
watch(
    () => form.ptm_localization,
    (on) => {
        if (on && !form.intensity.includes("PTM")) {
            form.intensity = "Prosit_2024_intensity_PTMs_gl";
            form.irt = "Prosit_2024_irt_PTMs_gl";
        }
        if (on && form.quant_precursor === "flashlfq") form.quant_precursor = "search_engine";
    }
);

// new PTMs start with their default residues
watch(
    () => [...form.loc_ptms],
    (ids) => {
        for (const id of ids) {
            if (!form.loc_residues[id]) form.loc_residues[id] = [...(store.meta.localizationPtms[id]?.residues ?? [])];
        }
    }
);

watch(
    () => form.loc_neutral_loss,
    (on) => {
        if (!on) form.loc_decoy_neutral_loss = false;
    }
);

watch(
    () => form.unitMassTolerance,
    (unit) => {
        if (unit !== "ppm") form.massToleranceAuto = false;
    }
);


const isValid = computed(() => {
    const hasSearchResults = searchResultFiles.value.length > 0;
    const hasSpectra = spectraFiles.value.length > 0;

    const searchResultsDone = hasSearchResults && searchResultFiles.value.every((f) => f.status === "done");
    const spectraDone = hasSpectra && spectraFiles.value.every((f) => f.status === "done");
    const fastaDone = !form.quantification || (fastaFiles.value.length > 0 && fastaFiles.value.every((f) => f.status === "done"));
    const ceValid = form.ce_range[0] < form.ce_range[1];

    return searchResultsDone && spectraDone && fastaDone && ceValid && !quantProblem.value && !localizationProblem.value;
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
        tag: form.tag,
        fdr_estimation_method: form.fdr_estimation_method,
        regressionMethod: form.regressionMethod,
        add_feature_cols: form.add_feature_cols,
        massTolerance: form.massToleranceAuto ? "auto" : form.massTolerance,
        unitMassTolerance: form.unitMassTolerance,
        ion_types: [...form.ion_types_list].sort().reverse().join(""),
        ce_alignment_options: {
            ce_range: form.ce_range,
            use_ransac_model: form.use_ransac_model,
        },
        quantification: form.quantification,
        percolatorSubsetMaxTrain: form.percolatorSubsetEnabled ? form.percolatorSubsetMaxTrain : null,
        chimericPeakSubtraction: form.chimericPeakSubtraction,
        ptm_localization: form.ptm_localization,
    };
    if (form.ptm_fdr_control !== "auto") {
        cfg["ptm_fdr_control"] = form.ptm_fdr_control === "on";
    }
    if (form.quantification) {
        const quant: Record<string, unknown> = { precursor: form.quant_precursor };
        if (isFlashlfq.value) {
            Object.assign(quant, {
                protein: form.quant_protein,
                normalization: form.quant_normalization,
                psmFdr: form.quant_psmFdr,
                peptideFdr: form.quant_peptideFdr,
                proteinFdr: form.quant_proteinFdr,
                ppmTolerance: form.quant_ppmTolerance,
                ppmFilter: form.quant_ppmFilterMode === "value" ? form.quant_ppmFilter : form.quant_ppmFilterMode,
            });
            if (form.quant_mbr) quant["mbrFdr"] = form.quant_mbrFdr;
        }
        cfg["quantificationOptions"] = quant;
        cfg["fastaDigestOptions"] = {
            digestion: form.digest_digestion,
            missedCleavages: form.digest_missedCleavages,
            minLength: form.digest_minLength,
            maxLength: form.digest_maxLength,
            enzyme: form.digest_enzyme,
            specialAas: form.digest_specialAas,
            db: "concat",
            fragmentation: "HCD",
        };
    }
    if (form.ptm_localization) {
        const loc: Record<string, unknown> = {
            ptms: Object.fromEntries(form.loc_ptms.map((id) => [id, form.loc_residues[id] ?? []])),
            neutral_loss: form.loc_neutral_loss,
            decoy_neutral_loss: form.loc_decoy_neutral_loss,
            decoy_selection: form.loc_decoy_selection,
            flr_neighbor_split: form.loc_flr_neighbor_split,
            loc_probability: form.loc_probability,
            decoy_weighting: form.loc_decoy_weighting,
            flr_groups: form.loc_flr_groups,
            flr_estimator: form.loc_flr_estimator,
            sdi_features: form.loc_sdi_features,
            ambiguity_check: form.loc_ambiguity_check,
        };
        if (decoyResidueList.value.length) loc["decoy_residues"] = decoyResidueList.value;
        if (form.loc_decoy_proxy) loc["decoy_proxy"] = form.loc_decoy_proxy;
        cfg["ptmLocalizationOptions"] = loc;
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
                            @upload="store.uploadRole('Rescoring', 'search_results')"
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
                            @upload="store.uploadRole('Rescoring', 'spectra')"
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
                            :items="['percolator']"
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
                            :disabled="form.massToleranceAuto"
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
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.tag"
                            :items="store.meta.tags"
                            label="Isobaric tag"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.instrument_type"
                            :items="['QE', 'LUMOS', 'TIMSTOF', 'SCIEXTOF']"
                            label="Instrument type"
                            density="compact"
                            variant="outlined"
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
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.massToleranceAuto"
                            label="Auto mass tolerance (per file, ppm)"
                            color="primary"
                            hide-details
                            :disabled="form.unitMassTolerance !== 'ppm'"
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.chimericPeakSubtraction"
                            label="Chimeric peak subtraction"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.percolatorSubsetEnabled"
                            label="Train Percolator on a subset"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-text-field
                            v-model.number="form.percolatorSubsetMaxTrain"
                            label="Max. PSMs for Percolator training"
                            type="number"
                            min="1"
                            density="compact"
                            variant="outlined"
                            :disabled="!form.percolatorSubsetEnabled"
                            hint="Speeds up large inputs; all PSMs are still scored"
                            persistent-hint
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.ptm_fdr_control"
                            :items="[
                                { title: 'Auto (on for PTM models without localization)', value: 'auto' },
                                { title: 'On', value: 'on' },
                                { title: 'Off', value: 'off' },
                            ]"
                            label="Separate FDR for modified PSMs"
                            density="compact"
                            variant="outlined"
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
                <v-row>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.quant_precursor"
                            :items="precursorItems"
                            label="Precursor quantification"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="text-body-2 text-medium-emphasis d-flex align-center">
                        <span v-if="isFlashlfq">
                            Measures the MS1 signal of the rescored PSMs with FlashLFQ, protein groups by
                            picked-group-FDR.
                        </span>
                        <span v-else>
                            Updates the search engine's quantification with the rescored PSMs (picked-group-FDR).
                            Upload msms.txt + evidence.txt (MaxQuant) or results.sage.tsv + lfq.tsv (Sage).
                        </span>
                    </v-col>
                </v-row>
                <v-row v-if="isFlashlfq">
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.quant_protein"
                            :items="store.meta.quantification.protein"
                            label="Protein intensities"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.quant_normalization"
                            :items="store.meta.quantification.normalization"
                            label="Normalization"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.quant_ppmTolerance"
                            label="MS1 tolerance (ppm)"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.quant_psmFdr"
                            label="PSM q-value cutoff"
                            type="number"
                            step="0.01"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.quant_peptideFdr"
                            label="Peptide q-value cutoff"
                            type="number"
                            step="0.01"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.quant_proteinFdr"
                            label="Protein group q-value cutoff"
                            type="number"
                            step="0.01"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.quant_ppmFilterMode"
                            :items="[
                                { title: 'Auto (90th percentile per raw file)', value: 'auto' },
                                { title: 'Fixed window (ppm)', value: 'value' },
                                { title: 'Off (keep all peaks)', value: 'none' },
                            ]"
                            label="Peak mass error filter"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-text-field
                            v-model.number="form.quant_ppmFilter"
                            label="Mass error window (ppm)"
                            type="number"
                            density="compact"
                            variant="outlined"
                            :disabled="form.quant_ppmFilterMode !== 'value'"
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.quant_mbr"
                            label="Match between runs"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-text-field
                            v-model.number="form.quant_mbrFdr"
                            label="MBR q-value cutoff (PIP)"
                            type="number"
                            step="0.01"
                            density="compact"
                            variant="outlined"
                            :disabled="!form.quant_mbr"
                        />
                    </v-col>
                </v-row>
                <div class="text-body-2 mt-2 mb-2 font-weight-medium">
                    In-silico digestion (as in the search)
                </div>
                <v-row>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.digest_enzyme"
                            :items="store.meta.enzymes"
                            label="Enzyme"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-select
                            v-model="form.digest_digestion"
                            :items="['full', 'semi', 'none']"
                            label="Digestion"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model="form.digest_specialAas"
                            label="Special AAs (cleavage sites)"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.digest_missedCleavages"
                            label="Missed cleavages"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.digest_minLength"
                            label="Min. peptide length"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                    <v-col cols="12" md="4">
                        <v-text-field
                            v-model.number="form.digest_maxLength"
                            label="Max. peptide length"
                            type="number"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                </v-row>
                <div class="text-body-2 mb-2 font-weight-medium">
                    FASTA file (for quant) *
                </div>
                <FileUpload
                    role="fasta"
                    :accept="['.fasta', '.fa', '.faa']"
                    :multiple="false"
                    @add="onFilesAdded('fasta', $event)"
                    @upload="store.uploadRole('Rescoring', 'fasta')"
                />
                <v-alert v-if="quantProblem" type="warning" variant="tonal" density="compact" class="mt-4">
                    {{ quantProblem }}
                </v-alert>
            </v-card-text>
        </v-card>

        <!-- PTM localization -->
        <v-card variant="outlined" class="mb-4">
            <v-card-title
                class="d-flex align-center text-body-1 font-weight-medium"
            >
                <span>PTM Localization</span>
                <v-switch
                    v-model="form.ptm_localization"
                    color="primary"
                    hide-details
                    class="ml-4"
                />
            </v-card-title>
            <v-card-text v-if="form.ptm_localization">
                <div class="d-flex align-center text-body-2 text-medium-emphasis mb-4">
                    <span>
                        Scores all site localizations of the selected PTMs and estimates the false localization rate
                        (FLR) with decoy sites. Needs a PTM intensity model. Defaults are the settings of the
                        published phospho benchmarks.
                    </span>
                    <v-btn variant="text" size="small" color="primary" class="ml-2" @click="resetLocalizationDefaults">
                        Reset to defaults
                    </v-btn>
                </div>
                <v-row>
                    <v-col cols="12">
                        <v-select
                            v-model="form.loc_ptms"
                            :items="ptmItems"
                            label="PTMs to localize"
                            density="compact"
                            variant="outlined"
                            multiple
                            chips
                        />
                    </v-col>
                    <v-col v-for="id in form.loc_ptms" :key="id" cols="12" md="6">
                        <v-select
                            v-model="form.loc_residues[id]"
                            :items="AMINO_ACIDS"
                            :label="`Possible residues of ${store.meta.localizationPtms[id]?.name ?? 'UNIMOD:' + id}`"
                            density="compact"
                            variant="outlined"
                            multiple
                            chips
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.loc_neutral_loss"
                            label="Neutral loss feature"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.loc_sdi_features"
                            label="Site-determining ion features"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6" class="d-flex align-center">
                        <v-switch
                            v-model="form.loc_ambiguity_check"
                            label="Rank ambiguous localizations last"
                            color="primary"
                            hide-details
                        />
                    </v-col>
                    <v-col cols="12" md="6">
                        <v-select
                            v-model="form.loc_flr_estimator"
                            :items="[
                                { title: 'Hybrid (decoy + model, default)', value: 'hybrid' },
                                { title: 'Decoy only', value: 'decoy' },
                            ]"
                            label="FLR estimator"
                            density="compact"
                            variant="outlined"
                        />
                    </v-col>
                </v-row>
                <v-expansion-panels variant="accordion" class="mt-2">
                    <v-expansion-panel title="Advanced localization options">
                        <v-expansion-panel-text>
                            <v-row>
                                <v-col cols="12" md="6">
                                    <v-text-field
                                        v-model="form.loc_decoy_residues"
                                        label="Decoy residues (by priority)"
                                        placeholder="empty: no decoys (no FLR estimate)"
                                        density="compact"
                                        variant="outlined"
                                        hint="Residues that cannot carry the PTMs, comma separated"
                                        persistent-hint
                                    />
                                </v-col>
                                <v-col cols="12" md="6">
                                    <v-select
                                        v-model="form.loc_decoy_selection"
                                        :items="[
                                            { title: 'Residue (one decoy residue per peptide, default)', value: 'residue' },
                                            { title: 'Adjacent (also near possible sites, experimental)', value: 'adjacent' },
                                        ]"
                                        label="Decoy selection"
                                        density="compact"
                                        variant="outlined"
                                    />
                                </v-col>
                                <v-col cols="12" md="6">
                                    <v-select
                                        v-model="form.loc_decoy_proxy"
                                        :items="[
                                            { title: `Default (${possibleSites[0] ?? 'first possible site'})`, value: '' },
                                            { title: 'Search engine frequencies', value: 'search' },
                                            ...possibleSites.map((aa) => ({ title: aa, value: aa })),
                                        ]"
                                        label="Decoy prediction proxy"
                                        density="compact"
                                        variant="outlined"
                                    />
                                </v-col>
                                <v-col cols="12" md="6" class="d-flex align-center">
                                    <v-switch
                                        v-model="form.loc_decoy_neutral_loss"
                                        label="Neutral loss on decoy residues"
                                        color="primary"
                                        hide-details
                                        :disabled="!form.loc_neutral_loss"
                                    />
                                </v-col>
                                <v-col cols="12" md="6">
                                    <v-select
                                        v-model="form.loc_probability"
                                        :items="[
                                            { title: 'PEP (default)', value: 'pep' },
                                            { title: 'Score (robust to collapsed PEPs)', value: 'score' },
                                        ]"
                                        label="Localization probability"
                                        density="compact"
                                        variant="outlined"
                                    />
                                </v-col>
                                <v-col cols="12" md="6">
                                    <v-select
                                        v-model="form.loc_decoy_weighting"
                                        :items="[
                                            { title: 'Count (default)', value: 'count' },
                                            { title: 'Distance', value: 'distance' },
                                        ]"
                                        label="Decoy weighting"
                                        density="compact"
                                        variant="outlined"
                                    />
                                </v-col>
                                <v-col cols="12" md="6">
                                    <v-select
                                        v-model="form.loc_flr_groups"
                                        :items="[
                                            { title: 'One estimate (default)', value: 'none' },
                                            { title: 'Per residue', value: 'residue' },
                                        ]"
                                        label="FLR groups"
                                        density="compact"
                                        variant="outlined"
                                    />
                                </v-col>
                                <v-col cols="12" md="6" class="d-flex align-center">
                                    <v-switch
                                        v-model="form.loc_flr_neighbor_split"
                                        label="Separate FLR next to possible sites"
                                        color="primary"
                                        hide-details
                                    />
                                </v-col>
                            </v-row>
                        </v-expansion-panel-text>
                    </v-expansion-panel>
                </v-expansion-panels>
                <v-alert v-if="localizationProblem" type="warning" variant="tonal" density="compact" class="mt-4">
                    {{ localizationProblem }}
                </v-alert>
            </v-card-text>
        </v-card>


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
