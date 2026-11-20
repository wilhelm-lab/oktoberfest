<script setup lang="ts">
import { ref, onMounted, watch } from "vue";
import { useRoute } from "vue-router";
import UploadProgressOverlay from "@/components/UploadProgressOverlay.vue";
import { getSession } from "@/api/client";

// show when this browser is signed in as admin (hosted mode); re-checked after leaving the admin page
const route = useRoute();
const isAdmin = ref(false);

async function refreshSession() {
    try {
        isAdmin.value = (await getSession()).is_admin;
    } catch {
        isAdmin.value = false;
    }
}

onMounted(refreshSession);
watch(() => route.name, (name, previous) => {
    if (previous === "admin") refreshSession();
});
</script>

<template>
    <v-app>
        <v-app-bar color="primary" elevation="2">
            <v-app-bar-title>
                <router-link
                    to="/"
                    style="
                        color: white;
                        text-decoration: none;
                        font-weight: 700;
                    "
                >
                    Oktoberfest
                </router-link>
            </v-app-bar-title>
            <template #append>
                <v-chip
                    v-if="isAdmin"
                    :to="{ name: 'admin' }"
                    color="white"
                    variant="outlined"
                    prepend-icon="mdi-shield-account"
                    class="mr-2"
                    >Admin</v-chip
                >
                <v-btn
                    variant="text"
                    color="white"
                    href="https://github.com/wilhelm-lab/oktoberfest"
                    target="_blank"
                    prepend-icon="mdi-github"
                    >Repo</v-btn
                >
                <v-btn
                    variant="text"
                    color="white"
                    href="https://oktoberfest.readthedocs.io"
                    target="_blank"
                    prepend-icon="mdi-book-open-variant"
                    >Docs</v-btn
                >
            </template>
        </v-app-bar>

        <v-main>
            <router-view />
        </v-main>

        <v-footer color="primary" class="text-caption text-center py-2">
            <span style="color: rgba(255, 255, 255, 0.7)">
                Oktoberfest Web UI ·
                <a
                    href="https://github.com/wilhelm-lab/oktoberfest"
                    target="_blank"
                    style="color: #bfda79"
                    >GitHub</a
                >
                ·
                <a
                    href="https://oktoberfest.readthedocs.io"
                    target="_blank"
                    style="color: #bfda79"
                    >Docs</a
                >
            </span>
        </v-footer>
        <UploadProgressOverlay />
    </v-app>
</template>
