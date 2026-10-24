import { describe, it, expect } from "vitest";
import { mount } from "@vue/test-utils";
import { createVuetify } from "vuetify";
import * as components from "vuetify/components";
import * as directives from "vuetify/directives";
import ConfigSummary from "../ConfigSummary.vue";

const vuetify = createVuetify({ components, directives });

describe("ConfigSummary", () => {
    it("renders the JSON of the config", () => {
        const config = {
            type: "SpectralLibraryGeneration",
            models: {
                intensity: "Prosit_2020_intensity_HCD",
                irt: "Prosit_2019_irt",
            },
        };
        const wrapper = mount(ConfigSummary, {
            props: { config },
            global: { plugins: [vuetify] },
        });
        expect(wrapper.text()).toContain("SpectralLibraryGeneration");
        expect(wrapper.text()).toContain("Prosit_2020_intensity_HCD");
    });

    it("pretty-prints the JSON", () => {
        const config = { type: "Rescoring", numThreads: 4 };
        const wrapper = mount(ConfigSummary, {
            props: { config },
            global: { plugins: [vuetify] },
        });
        expect(wrapper.html()).toContain('"numThreads": 4');
    });
});
