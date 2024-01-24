<template>

  <v-flex xs9>
    <v-card class="grey lighten-4">
      <v-toolbar color="primary" dark flat dense>
        <v-toolbar-title>Change Log</v-toolbar-title>
      </v-toolbar>
      <v-card-text>



        <v-row justify="center">
                        <v-expansion-panel accordion :value="0">
                          <v-expansion-panel-content v-for="post in logs"
                                                     :key="post.title">
                            <div class="posts" slot="header"><b>{{ post.title }}</b></div>
                               <v-card>
                                <v-card-text>
                                  &nbsp;&nbsp;&nbsp;&nbsp;{{ post.date }}<br>
                                  &nbsp;&nbsp;&nbsp;&nbsp;{{ post.text }}
                                </v-card-text>
                               </v-card>
                          </v-expansion-panel-content>

                        </v-expansion-panel>

        </v-row>

      </v-card-text>

    </v-card>
  </v-flex>
</template>


<script>

import axios from "axios";

export default {
  name: 'LOG',
  methods: {
    newTask: function(){
      this.$router.push('/prosit/');
    }
  },
  mounted() {
    let self = this;
    //self.logs = this.$store.state.posts
    let logsAPI = '/prosit/api/changelog.xsjs';
    axios.get(logsAPI).then(function(response){
      self.logs = response.data.messages;
    });
  },
  data: () => ({
    logs: []
  }),
}


</script>
<style scoped>
.posts {
  list-style: none;
  text-align: left;
}
</style>
