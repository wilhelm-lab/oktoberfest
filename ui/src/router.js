import Vue from 'vue'
import Router from 'vue-router'
import TaskSelection from './views/TaskSelection.vue'
import TaskStatus from './views/TaskStatus.vue'
import FAQ from './views/FAQ.vue'
import Libraries from './views/Libraries.vue'
import Changelog from "@/views/Changelog";



Vue.use(Router)

export default new Router({
  routes: [
    {
      path: '/prosit/',
      name: 'home',
      component: TaskSelection
    },
    {
      path: '/prosit/task/:taskid',
      name: 'task',
      props: true,
      component: TaskStatus
    },
    {
      path: '/prosit/faq/',
      name: 'faq',
      component: FAQ
    },
    {
      path: '/prosit/log/',
      name: 'log',
      component: Changelog
    },
    {
      path: '/prosit/libraries/',
      name: 'libraries',
      component: Libraries
    }
  ],
  mode: 'history'
})
