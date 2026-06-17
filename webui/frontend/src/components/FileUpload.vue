<script setup lang="ts">
import { ref } from "vue";

const props = defineProps<{
    role: string;
    accept: string[];
    multiple?: boolean;
    sizeLimit?: number; // bytes
}>();

const emit = defineEmits<{
    (e: "add", files: File[]): void;
    (e: "remove", index: number): void;
}>();

const fileInput = ref<HTMLInputElement | null>(null);
const dragover = ref(false);
const selectedFiles = ref<{ file: File; error?: string }[]>([]);

const acceptStr = props.accept.join(",");

function validate(file: File): string | null {
    const ext = "." + file.name.split(".").pop()!.toLowerCase();
    if (
        props.accept.length &&
        !props.accept.map((a) => a.toLowerCase()).includes(ext)
    ) {
        return `Extension ${ext} not allowed. Allowed: ${props.accept.join(
            ", "
        )}`;
    }
    if (props.sizeLimit && file.size > props.sizeLimit) {
        return `File too large (max ${(props.sizeLimit / 1e9).toFixed(1)} GB)`;
    }
    return null;
}

function addFiles(files: FileList | File[]) {
    const arr = Array.from(files);
    const newEntries = arr.map((f) => ({
        file: f,
        error: validate(f) ?? undefined,
    }));
    if (!props.multiple) selectedFiles.value = [];
    selectedFiles.value.push(...newEntries);
    const valid = newEntries.filter((e) => !e.error).map((e) => e.file);
    if (valid.length) emit("add", valid);
}

function onInputChange(e: Event) {
    const t = e.target as HTMLInputElement;
    if (t.files) addFiles(t.files);
    t.value = "";
}

function onDrop(e: DragEvent) {
    dragover.value = false;
    if (e.dataTransfer?.files) addFiles(e.dataTransfer.files);
}

function remove(i: number) {
    selectedFiles.value.splice(i, 1);
    emit("remove", i);
}

function formatSize(bytes: number): string {
    if (bytes < 1024) return bytes + " B";
    if (bytes < 1024 ** 2) return (bytes / 1024).toFixed(1) + " KB";
    if (bytes < 1024 ** 3) return (bytes / 1024 ** 2).toFixed(1) + " MB";
    return (bytes / 1024 ** 3).toFixed(2) + " GB";
}
</script>

<template>
    <div>
        <v-card
            :class="['drop-zone pa-4 text-center', { dragover: dragover }]"
            variant="outlined"
            @dragover.prevent="dragover = true"
            @dragleave.prevent="dragover = false"
            @drop.prevent="onDrop"
            @click="fileInput?.click()"
            style="cursor: pointer; border: 2px dashed rgba(14, 45, 78, 0.3)"
        >
            <v-icon size="40" color="primary">mdi-cloud-upload</v-icon>
            <div class="mt-2 text-body-2">
                Drop files here or <strong>click to browse</strong>
            </div>
            <div class="text-caption text-grey">
                Accepted: {{ accept.join(", ") || "any" }}
            </div>
            <input
                ref="fileInput"
                type="file"
                :accept="acceptStr"
                :multiple="multiple"
                style="display: none"
                @change="onInputChange"
            />
        </v-card>

        <v-list v-if="selectedFiles.length" density="compact" class="mt-2">
            <v-list-item
                v-for="(entry, i) in selectedFiles"
                :key="i"
                :subtitle="entry.error ?? formatSize(entry.file.size)"
                :class="entry.error ? 'text-error' : ''"
            >
                <template #prepend>
                    <v-icon :color="entry.error ? 'error' : 'success'">
                        {{
                            entry.error ? "mdi-alert-circle" : "mdi-file-check"
                        }}
                    </v-icon>
                </template>
                <template #title>{{ entry.file.name }}</template>
                <template #append>
                    <v-btn
                        icon="mdi-close"
                        size="small"
                        variant="text"
                        @click.stop="remove(i)"
                    />
                </template>
            </v-list-item>
        </v-list>
    </div>
</template>

<style scoped>
.drop-zone.dragover {
    border-color: #0e2d4e !important;
    background: rgba(191, 218, 121, 0.1);
}
</style>
