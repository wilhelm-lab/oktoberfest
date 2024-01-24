<template>
  <v-app>
  <v-container fluid grid-list-lg>
    <v-layout row wrap class="justify-center">
      
      <v-flex xs9>
        <v-alert :value='underMaintenance' type="error" color="#D77373">
          Prosit is currently under maintenance. Please excuse the inconvenience.
        </v-alert>
      </v-flex>

      <v-flex xs9>
        <v-alert type="warning" color="#d6983c" value="true">
          We now offer two new Prosit TMT models that will soon be published. One is for fragment intensities prediction <b>(Prosit_TMT_intensity_2021)</b>
          and the other is for iRT prediction <b>(Prosit_TMT_irt_2021)</b>. The intensity model works with both <b>CID</b> and <b>HCD</b> fragmentation methods but you need
          to add fragmentation column to the input. We assume all the sequences are fully labeled and you don't need to add the
          tmt modification explicitly in your input files.
        </v-alert>
      </v-flex>

      <v-flex xs9>

        <v-card>
          
          <v-toolbar color="primary" dark flat>
            <v-toolbar-title @click="showHome" :style="{ cursor: 'pointer'}">  
              Prosit
            </v-toolbar-title>
            <img src="@/assets/prosit_pixel.png" alt="" style="width:35px;height:50px;margin-left:20px;">
            <!--<v-chip class="ma-2" color="primary" label text-color="white">No. of jobs waiting in queue: {{runJobs}}</v-chip>-->
            <v-spacer></v-spacer>
            <!--<v-btn  outline color="white" @click="showLOG">Changelog</v-btn>-->
            <v-btn  outline :disabled='underMaintenance' color="white" @click="showHome">predict</v-btn>
            <v-btn  outline color="white" @click="showLibraries">libraries</v-btn>
            <v-btn  outline color="white" @click="showFAQ">faq</v-btn>
            <v-btn  outline :disabled='underMaintenance' color="white" @click="toggleTaskDialog">status</v-btn>

            <v-btn icon href="http://github.com/kusterlab/prosit" target="_blank">
              <v-icon>fab fa-github</v-icon>
            </v-btn>
            <v-btn icon href="mailto: prosit.proteomics.wzw@tum.de" target="_blank">
              <v-icon>fas fa-question-circle</v-icon>
            </v-btn>
          </v-toolbar>

          <v-card-text id='abstract' class="grey lighten-4">
            Prosit offers high quality MS2 predicted spectra for any organism and protease as well as iRT prediction.
            Prosit is part of the ProteomeTools (<a href='http://www.proteometools.org' target="_blank">www.proteometools.org/</a>)
            project and was trained on the project's high quality synthetic dataset.
            When using Prosit is helpful for your research, please cite "Gessulat, Schmidt et al. 2019" <a href='https://doi.org/10.1038/s41592-019-0426-7' target="_blank">DOI 10.1038/s41592-019-0426-7</a>.
          </v-card-text>

        </v-card>
      </v-flex>

      <router-view :hidden='hideRouter'/> 

    </v-layout>
  </v-container>

  <v-dialog v-model="displayTaskDialog" max-width="300px">
  <v-card>
    <v-card-text>
      <v-text-field label="Task ID" v-model="taskid" clearable @keyup.enter="showTaskStatus"></v-text-field>
    </v-card-text>
    <v-card-actions>
      <v-spacer></v-spacer>
      <v-btn :disabled="!taskid" small fab color="primary" @click.native="showTaskStatus"><v-icon>done</v-icon></v-btn>
    </v-card-actions>
  </v-card>
</v-dialog>

  </v-app>
</template>

<script>
import router from '@/router'
import axios from 'axios'



export default {
  name: 'app',
  data: () => ({
    displayTaskDialog: false,
    taskid: null,
    underMaintenance: false,
    runJobs: -1,
  }),
  methods: {
    toggleTaskDialog: function(){
      this.displayTaskDialog = ! this.displayTaskDialog;
    },
    showHome: function(){
      router.push('/prosit/')
    },
    showLOG: function(){
      router.push('/prosit/log/')
    },
    showFAQ: function(){
      router.push('/prosit/faq/')
    },
    showLibraries: function(){
      router.push('/prosit/libraries/')
    },
    showTaskStatus: function(){
      this.toggleTaskDialog();
      router.push('/prosit/task/' + this.taskid)
    },
  },
  computed: {
    hidePredict: function () {
      return this.$route.name == 'home' || this.underMaintenance;
    },
    hideRouter: function () {
      if(this.underMaintenance){
        return this.$route.name != 'faq'
      }
      return false;
    }
  },
  mounted() {
    let statusUrl = '/prosit/api/backend.xsjs';
    let self = this;
    // expects an json object: {"underMaintenance": true} or false.
    axios.get(statusUrl).then(function (response){   
        self.underMaintenance = response.data.underMaintenance == "true";
    });

    let runJobsUrl = '/prosit/api/get_running_jobs.xsjs';
    axios.get(runJobsUrl).then(function (response){
      self.runJobs = response.data.number_of_jobs;
    });
  }
}
</script>

<style>
</style>