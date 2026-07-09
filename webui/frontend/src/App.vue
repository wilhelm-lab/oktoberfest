<script setup lang="ts">
import UploadProgressOverlay from "@/components/UploadProgressOverlay.vue";
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
