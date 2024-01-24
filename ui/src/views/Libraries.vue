<template>
  <v-layout class="justify-center">

    <v-flex xs9>
      <v-card >
        <v-toolbar color="primary" dark flat dense>
          <v-toolbar-title>Spectral Libraries</v-toolbar-title>
        </v-toolbar>
        <v-card-text>
          <p>
          The following are ready-made predicted spectral library to use with Encyclopedia for a new DIA-only spectral library generation workflow. <br />
          For more information see "Searle et. al. 2019" <a href="https://doi.org/10.1101/682245 ">DOI: 10.1101/682245</a> <br />
          <b>NOTE:</b> some files are exceedingly large. Initiating the download may take a while.
          </p>
          <ul>
            <li v-for="item in librariesInfo" v-bind:key="item.Header"><b>{{item.Header}}</b><br />
              <div v-for="item in item.File"  v-bind:key="item.Path">
                <a  :href="getUrl(item.Path)">{{item.Description}}</a> <br />
              </div>
            </li>
          </ul>
        </v-card-text>
      </v-card>
    </v-flex>
  </v-layout>
</template>


<script>
import axios from 'axios'

export default {
  name: 'Libraries',
  mounted() {
    let self = this;
    let librariesAPI = '/prosit/api/libraries.xsjs';
    axios.get(librariesAPI).then(function(response){
      self.librariesInfo = response.data.Spectral_library;
    });
  },
  data: () => ({
    librariesInfo: []
  }),
  methods: {
    getUrl: function(fileName) {
      return "/prosit/api/libraries.xsjs?file=" + fileName;
    }
  }
}
</script>