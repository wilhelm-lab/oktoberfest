<template>
        <v-container>

        <v-snackbar v-model="fileIsTooBig" :timeout=3000 top>
        The file exceeds the size limit of {{sizelimit}}MB.
        </v-snackbar>

        <v-text-field 
            readonly 
            loading
            :hint='hinttext'
            persistent-hint
            @click="selectFile"
            @click:clear="resetUpload"
            v-model="filename" 
            append-outer-icon="clear"
            @click:append-outer="clear"
            class='ma-0'> 
          <v-tooltip top slot='prepend'>
            <v-btn :disabled='progress != 0' slot='activator' fab small color='primary' @click.native="selectFile"> 
              <v-icon>cloud_upload</v-icon>
            </v-btn>
            {{this.filetype}}
          </v-tooltip>
          <v-progress-linear
                color="accent"
                height="3"
                slot="progress"
                :value="progress"
          ></v-progress-linear>
          
        </v-text-field>
        <!-- <input hidden type="file" ref="file" v-on:change="uploadFile()" :accept='filesuffix'/> -->
        <input hidden type="file" ref="file" v-on:change="verifyUpload()" :accept='filesuffix'/>
        </v-container>


</template>

<script>
import axios from 'axios'
import store from '@/store'


export default {
  props: {
    'filetype': String, 
    'taskid': String, 
    'hinttext': String, 
    'filesuffix': String,
  },
  data () {
    return {
      loading: false,
      progress: 0,
      baseurl: '/prosit/api/upload.xsjs',
      file: null,
      filename: this.filetype,
      fileIsTooBig: false,
      requesttoken: axios.CancelToken.source()
    }
  },

  methods: {
    clear(){
      this.requesttoken.cancel();
      this.resetUpload();

    },
    resetUpload(){
      this.progress = 0;
      this.loading = false;
      this.filename = this.filetype;
      this.$refs.file.value = "";
      store.commit('reportUpload', {fileType: this.filetype, isSuccessful: false});
    },
    selectFile(){
      this.$refs.file.click()
    },
    verifyUpload(){
      var file = this.$refs.file.files[0];
      var filesizeMB = (file.size/1024)/1024;
      this.fileIsTooBig = this.sizelimit > -1 && this.sizelimit < filesizeMB;
      if(this.fileIsTooBig){
        this.clear();
      }else{
        this.uploadFile();
      }
    },
    uploadFile(){
      this.requesttoken.cancel();
      this.progress = 0;

      this.file = this.$refs.file.files[0];
      this.filename = this.file.name;
      let formData = new FormData();
      formData.set('type', this.filetype);
      formData.append('file', this.file);
      this.requesttoken = axios.CancelToken.source();
      let self = this;

      axios.post(
        this.baseurl,
        formData,
        {
          params: { datasetId: this.taskid },
          headers: { 'Content-Type': 'multipart/form-data' },
          onUploadProgress: function( progressEvent ) {
            this.progress = parseInt( Math.round( ( progressEvent.loaded * 100 ) / progressEvent.total ) );
          }.bind(this),
          cancelToken: this.requesttoken.token
        }
      ).then(function(){ 
        // sucess
        store.commit('reportUpload', {fileType: self.filetype, isSuccessful: true});
      }).catch(
        // error
        self.resetUpload 
      );
    },
  }
}
</script>


