<template>
  <v-layout class="justify-center">

    <v-flex xs9>
      <v-card class="grey lighten-4">
        <v-toolbar color="primary" dark flat dense>
          <v-toolbar-title>Task {{taskid}}</v-toolbar-title>
        </v-toolbar>
        <v-card-text id='abstract' >
          {{statusText}}


          <div :hidden="prosit_stop_code <= 0">
          <v-container mb-3 >
          <v-expansion-panel>
          <v-expansion-panel-content>
            <div slot="header">Error log</div>
            <v-card>
              <v-card-text>{{prosit_error_message}}</v-card-text>
            </v-card>
          </v-expansion-panel-content>
          </v-expansion-panel>
          </v-container>
          </div>


        </v-card-text>

        <div :hidden="prosit_stop_code !== 0">
        <v-container >
            <v-btn small color="primary"
                  :href="downloadUrl" target="_blank">download</v-btn>
        </v-container>
        </div>

        <div >


        </div>

      </v-card>
    </v-flex>
  </v-layout>
</template>


<script>
import axios from 'axios'


export default {
  name: 'TaskStatus',
  components: {},
  props: ['taskid'],
  data: () => ({
    url: '/prosit/api/status.xsjs',
    download: "/prosit/api/download.xsjs?datasetId=",
    prosit_stop_code: -2,
    prosit_error_message: ""
  }),
  methods: {
    handleStatus: function(response){
      this.prosit_stop_code = parseInt(response.data.prosit_stop_code); 
      if(response.data.hasOwnProperty('prosit_error_message')){
        this.prosit_error_message = response.data.prosit_error_message.trim();
      }   
    }
  },
  computed: {
    downloadUrl: function () {
      return this.download + this.taskid;
    },
    statusText: function () {
      if(this.prosit_stop_code == -2) {
        return "This task ID is unkown. Please check if your task has another URL or re-submit your task."
      }
      else if(this.prosit_stop_code == -1) {
        return "This task is in progress. Tasks may take several hours for full proteomes depending on system load. " +
        "Please note down your Task ID or save this URL to check back later. " +
        "You can download the results here upon completion. " +
        "Resubmitting tasks will not lead to faster results. "
      }
      else if(this.prosit_stop_code == 0) {
        return "Your files are ready."
      }
      else {
        return "An error occured. Status code: " + this.prosit_stop_code
      }
    }
  },
  mounted() {
    axios.get(this.url, {params: {datasetId: this.taskid}}).then(this.handleStatus);
  }
}
</script>

