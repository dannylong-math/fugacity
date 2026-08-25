import { defineConfig, doxygen, markdown } from "sourcey";

export default defineConfig({
  name: "Fugacity",
  repo: "https://github.com/dannylong-math/fugacity",
  prettyUrls: false,
  theme: {
    colors: {
      primary: "#0F5A78",
      light: "#38BDF8",
      dark: "#0C4A6E",
    },
    fonts: {
      sans: "Inter",
    },
    layout: {
      content: "112rem",
    },
    css: ["./custom.css"],
  },
  navigation: {
    tabs: [
      {
        tab: "Guides",
        slug: "",
        source: markdown({
          groups: [
            {
              group: "Getting Started",
              pages: ["introduction", "getting-started", "tutorial"],
            },
            {
              group: "Extending Fugacity",
              pages: ["implementing-a-new-eos", "concepts"],
            },
          ],
        }),
      },
      {
        tab: "C++ API",
        slug: "api",
        source: doxygen({
          xml: "../build/doxygen/xml",
          language: "cpp",
          index: "flat",
        }),
      },
    ],
  },
});
