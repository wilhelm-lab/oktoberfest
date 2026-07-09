import { createRouter, createWebHistory } from "vue-router";

const router = createRouter({
    history: createWebHistory(import.meta.env.BASE_URL),
    routes: [
        {
            path: "/",
            name: "home",
            component: () => import("./views/HomeView.vue"),
        },
        {
            path: "/rescoring",
            name: "rescoring",
            component: () => import("./views/RescoringView.vue"),
        },
        {
            path: "/ce-calibration",
            name: "ce",
            component: () => import("./views/CeCalibrationView.vue"),
        },
        {
            path: "/spectral-library",
            name: "speclib",
            component: () => import("./views/SpecLibView.vue"),
        },
        {
            path: "/jobs",
            name: "jobs",
            component: () => import("./views/JobsListView.vue"),
        },
        {
            path: "/jobs/:jobId",
            name: "job",
            component: () => import("./views/JobStatusView.vue"),
        },
    ],
});

export default router;
