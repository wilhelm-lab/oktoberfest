<script setup lang="ts">
import { ref, onMounted } from "vue";
import { useRouter } from "vue-router";
import { getSession, adminSignIn, signOut, type SessionInfo } from "@/api/client";

const router = useRouter();
const session = ref<SessionInfo | null>(null);
const key = ref("");
const busy = ref(false);
const error = ref("");

onMounted(async () => {
    try {
        session.value = await getSession();
    } catch (e: unknown) {
        error.value = e instanceof Error ? e.message : String(e);
    }
});

async function signIn() {
    error.value = "";
    busy.value = true;
    try {
        session.value = await adminSignIn(key.value);
        key.value = "";
        router.push({ name: "jobs" });
    } catch (e: unknown) {
        error.value = e instanceof Error ? e.message : String(e);
    } finally {
        busy.value = false;
    }
}

async function leave() {
    error.value = "";
    busy.value = true;
    try {
        session.value = await signOut();
    } catch (e: unknown) {
        error.value = e instanceof Error ? e.message : String(e);
    } finally {
        busy.value = false;
    }
}
</script>

<template>
    <v-container max-width="560" class="py-6">
        <h2 class="text-h5 font-weight-bold text-primary mb-4">Admin</h2>

        <v-alert v-if="error" type="error" class="mb-4">{{ error }}</v-alert>

        <v-alert v-if="session && session.mode !== 'hosted'" type="info">
            This web UI runs in local mode: everyone already sees all jobs, no admin sign-in is needed.
        </v-alert>

        <v-card v-else-if="session?.is_admin" variant="outlined">
            <v-card-text>
                This browser is signed in as admin: the job list shows the jobs of all users, and you can cancel any
                of them.
            </v-card-text>
            <v-card-actions>
                <v-btn color="primary" :to="{ name: 'jobs' }">All jobs</v-btn>
                <v-btn variant="text" :loading="busy" @click="leave">Sign out</v-btn>
            </v-card-actions>
        </v-card>

        <v-card v-else-if="session" variant="outlined">
            <v-card-text>
                <p class="mb-4">
                    Sign in with an admin key (one of the <code>GLOBAL_USERS</code> values of the server) to see
                    and manage the jobs of all users in this browser.
                </p>
                <v-form @submit.prevent="signIn">
                    <v-text-field
                        v-model="key"
                        label="Admin key"
                        type="password"
                        autocomplete="current-password"
                        variant="outlined"
                        density="comfortable"
                    />
                    <v-btn type="submit" color="primary" :loading="busy" :disabled="!key">Sign in</v-btn>
                </v-form>
            </v-card-text>
        </v-card>

        <v-skeleton-loader v-else type="card" />
    </v-container>
</template>
