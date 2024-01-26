<template>
  <v-layout class="justify-center">

    <v-flex xs9>
      <v-card class="grey lighten-4">
        <v-toolbar color="primary" dark flat dense>
          <v-toolbar-title>Task {{taskid}}</v-toolbar-title>
        </v-toolbar>
        <v-card-text id='abstract' >
          {{statusText}}


          <div :hidden="status !== 'FAILED'">
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

        <div :hidden="status !== 'DONE'">
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
    url: process.env.VUE_APP_API_URL +'/api/v1/jobStatus',
    download: process.env.VUE_APP_API_URL +"/api/v1/downloadResults",
    status: false,
    prosit_error_message: ""
  }),
  methods: {
    handleStatus: function(response){
      this.status = response.data.status; 
      if(response.data.hasOwnProperty('errorMessage')){
        this.prosit_error_message = response.data.errorMessage.trim();
      }   
    }
  },
  computed: {
    downloadUrl: function () {
      return this.download + '?taskId=' + this.taskid;
    },
    statusText: function () {
      if(this.status == 'UNKNOWN') {
        return "This task ID is unkown. Please check if your task has another URL or re-submit your task."
      }
      else if(this.status == "PENDING" || this.status == "RUNNING") {
        return "This task is in progress. Tasks may take several hours for full proteomes depending on system load. " +
        "Please note down your Task ID or save this URL to check back later. " +
        "You can download the results here upon completion. " +
        "Resubmitting tasks will not lead to faster results. "
      }
      else if(this.status == "DONE") {
        return "Your files are ready."
      }
      else {
        return "An error occured. Status code: " + this.prosit_stop_code
      }
    }
  },
  mounted() {
    axios.get(this.url, {params: {taskId: this.taskid}}).then(this.handleStatus);
  }
}
</script>

