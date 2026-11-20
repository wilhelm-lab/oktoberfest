import { createApp } from "vue";
import { createPinia } from "pinia";
import { createVuetify } from "vuetify";
import * as components from "vuetify/components";
import * as directives from "vuetify/directives";
import "vuetify/styles";
import "@mdi/font/css/materialdesignicons.css";

import App from "./App.vue";
import router from "./router";

const vuetify = createVuetify({
    components,
    directives,
    theme: {
        defaultTheme: "light",
        themes: {
            light: {
                colors: {
                    primary: "#0E2D4E",
                    secondary: "#BFDA79",
                    accent: "#BFDA79",
                    surface: "#ffffff",
                },
            },
        },
    },
});

const app = createApp(App);
app.use(createPinia());
app.use(router);
app.use(vuetify);
app.mount("#app");
