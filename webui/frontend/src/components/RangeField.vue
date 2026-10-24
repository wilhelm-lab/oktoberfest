<script setup lang="ts">
const props = defineProps<{
    modelValue: [number, number];
    label: string;
    min?: number;
    max?: number;
}>();
const emit = defineEmits<{
    (e: "update:modelValue", v: [number, number]): void;
}>();

function updateMin(v: string) {
    emit("update:modelValue", [Number(v), props.modelValue[1]]);
}
function updateMax(v: string) {
    emit("update:modelValue", [props.modelValue[0], Number(v)]);
}
</script>

<template>
    <div>
        <div class="text-body-2 mb-1">{{ label }}</div>
        <div class="d-flex ga-2 align-center">
            <v-text-field
                :model-value="modelValue[0]"
                label="Min"
                type="number"
                density="compact"
                variant="outlined"
                hide-details
                style="max-width: 120px"
                @update:model-value="updateMin"
            />
            <span class="text-body-2">–</span>
            <v-text-field
                :model-value="modelValue[1]"
                label="Max"
                type="number"
                density="compact"
                variant="outlined"
                hide-details
                style="max-width: 120px"
                @update:model-value="updateMax"
            />
        </div>
        <div
            v-if="modelValue[0] >= modelValue[1]"
            class="text-error text-caption mt-1"
        >
            Min must be less than Max
        </div>
    </div>
</template>
