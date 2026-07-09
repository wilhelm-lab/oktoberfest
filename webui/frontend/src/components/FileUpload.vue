<script setup lang="ts">
import { ref, computed } from "vue";
import { useJobStore } from "@/stores/job";

const props = defineProps<{
    role: string;
    accept: string[];
    multiple?: boolean;
    sizeLimit?: number; // bytes
}>();

const emit = defineEmits<{
    (e: "add", files: File[]): void;
    (e: "upload"): void;
}>();

const store = useJobStore();
const fileInput = ref<HTMLInputElement | null>(null);
const dragover = ref(false);

const storeFiles = computed(() => store.uploads[props.role] || []);

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
    const valid = [];
    
    if (!props.multiple) {
        // Clear if not multiple. To do this cleanly, we need to clear store array
        store.uploads[props.role] = [];
    }

    for (const f of arr) {
        const error = validate(f);
        if (error) {
            // We can push to store with error status, or ignore
            // For now, let's let the store handle valid files
            alert(error); // simple error display
        } else {
            valid.push(f);
        }
    }
    
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
    store.uploads[props.role].splice(i, 1);
}

function formatSize(bytes: number): string {
    if (bytes < 1024) return bytes + " B";
    if (bytes < 1024 ** 2) return (bytes / 1024).toFixed(1) + " KB";
    if (bytes < 1024 ** 3) return (bytes / 1024 ** 2).toFixed(1) + " MB";
    return (bytes / 1024 ** 3).toFixed(2) + " GB";
}

const isUploadingThisRole = computed(() => storeFiles.value.some(f => f.status === "uploading"));
const isAllUploaded = computed(() => storeFiles.value.length > 0 && storeFiles.value.every(f => f.status === "done"));

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

        <v-list v-if="storeFiles.length" density="compact" class="mt-2">
            <v-list-item
                v-for="(entry, i) in storeFiles"
                :key="i"
                :class="entry.status === 'error' ? 'text-error' : ''"
            >
                <template #prepend>
                    <v-icon :color="entry.status === 'done' ? 'success' : (entry.status === 'error' ? 'error' : 'primary')">
                        {{
                            entry.status === 'done' ? 'mdi-check-circle' :
                            entry.status === 'error' ? 'mdi-alert-circle' :
                            entry.status === 'uploading' ? 'mdi-loading mdi-spin' :
                            'mdi-file-outline'
                        }}
                    </v-icon>
                </template>
                <template #title>{{ entry.file.name }}</template>
                <template #subtitle>
                    <div v-if="entry.status === 'uploading'">
                        <v-progress-linear :model-value="entry.progress" color="primary" height="4" class="mt-1" />
                    </div>
                    <div v-else>
                        {{ entry.error ?? formatSize(entry.file.size) }}
                    </div>
                </template>
                <template #append>
                    <v-btn
                        icon="mdi-close"
                        size="small"
                        variant="text"
                        :disabled="entry.status === 'uploading'"
                        @click.stop="remove(i)"
                    />
                </template>
            </v-list-item>
        </v-list>
        
        <div v-if="storeFiles.length" class="d-flex justify-end mt-2">
            <v-btn 
                size="small" 
                color="primary" 
                variant="tonal"
                prepend-icon="mdi-cloud-upload"
                @click="emit('upload')" 
                :disabled="isAllUploaded || isUploadingThisRole" 
                :loading="isUploadingThisRole"
            >
                {{ isAllUploaded ? 'Uploaded' : 'Upload' }}
            </v-btn>
        </div>
    </div>
</template>

<style scoped>
.drop-zone.dragover {
    border-color: #0e2d4e !important;
    background: rgba(191, 218, 121, 0.1);
}
</style>
