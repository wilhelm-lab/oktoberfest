import Vue from 'vue'
import App from '@/App.vue'
import router from '@/router'
import store from '@/store'

import Vuetify from 'vuetify'
import 'vuetify/dist/vuetify.min.css'
import 'material-design-icons-iconfont/dist/material-design-icons.css'
import '@fortawesome/fontawesome-free/css/all.css'


Vue.use(Vuetify, {
  theme: {
    primary: "#0E2D4E",
    secondary: "#BFDA79",
    accent: "#BFDA79"
  }
})


Vue.config.productionTip = false

new Vue({
  router,
  store,
  render: h => h(App)
}).$mount('#app')
