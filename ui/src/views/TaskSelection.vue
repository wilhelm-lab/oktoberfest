<template>
  <v-layout class="justify-center">

    <v-flex xs9>
      <v-card>
        <v-tabs grow centered dark color="primary" slider-color="accent" v-model="tabid" @click.native="changeTaskType">

          <v-tab ripple key="CEA">CE Calibration</v-tab>
          <v-tab-item><CEA :taskid="taskid"></CEA></v-tab-item>

          <v-tab ripple key="SLG">Spectral Library</v-tab>
          <v-tab-item><SLG :taskid="taskid"></SLG></v-tab-item>

          <!--<v-tab ripple key="PRP">Predict Properties</v-tab>
          <v-tab-item><PRP :taskid="taskid"></PRP></v-tab-item>-->

          <v-tab ripple key="MQR">Rescoring</v-tab>
          <v-tab-item><MQR :taskid="taskid" ></MQR></v-tab-item>

        </v-tabs>
      </v-card>
    </v-flex>
  </v-layout>
</template>


<script>
import axios from 'axios'
import SLG from '@/components/SLG.vue'
import CEA from '@/components/CEA.vue'
import MQR from '@/components/MQR.vue'
import PRP from '@/components/PRP.vue'
import store from '@/store'



export default {
  name: 'TaskSelection',
  components: {
    SLG, CEA, MQR,PRP
  },
  data: () => ({
    taskid: null,
    tabid: 0
  }),
  methods: {
    changeTaskType: function(){
      let taskType = store.state.taskTypes[this.tabid];
      store.commit('changeTask', taskType);
    },
  },
  mounted() {
  }
}
</script>
