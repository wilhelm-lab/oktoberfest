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
      path: '/',
      name: 'home',
      component: TaskSelection
    },
    {
      path: '/task/:taskid',
      name: 'task',
      props: true,
      component: TaskStatus
    },
    {
      path: '/faq/',
      name: 'faq',
      component: FAQ
    },
    {
      path: '/log/',
      name: 'log',
      component: Changelog
    },
    {
      path: '/libraries/',
      name: 'libraries',
      component: Libraries
    }
  ],
  mode: 'history'
})
