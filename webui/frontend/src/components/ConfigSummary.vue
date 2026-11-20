<script setup lang="ts">
import { computed } from "vue";

const props = defineProps<{
    config: Record<string, unknown>;
}>();

const json = computed(() => JSON.stringify(props.config, null, 2));

function copy() {
    navigator.clipboard.writeText(json.value);
}
</script>

<template>
    <v-card variant="outlined">
        <v-card-title class="d-flex align-center justify-space-between">
            <span class="text-body-1 font-weight-medium">Config preview</span>
            <v-btn
                size="small"
                variant="text"
                prepend-icon="mdi-content-copy"
                @click="copy"
            >
                Copy
            </v-btn>
        </v-card-title>
        <v-card-text>
            <pre
                class="text-caption"
                style="
                    overflow: auto;
                    max-height: 300px;
                    background: #f5f5f5;
                    padding: 8px;
                    border-radius: 4px;
                "
                >{{ json }}</pre
            >
        </v-card-text>
    </v-card>
</template>
