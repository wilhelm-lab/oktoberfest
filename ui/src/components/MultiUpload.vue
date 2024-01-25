<template>
  <v-container>
    <v-text-field 
            readonly 
            :hint='hinttext'
            persistent-hint
            @click="$refs.upload.$el.lastChild.click()"
            v-model="filetype" 
            class='ma-0'> 
      <v-tooltip top slot='prepend'>
        <v-btn slot='activator' fab small color='primary' @click.native="$refs.upload.$el.lastChild.click()">
          <v-icon>cloud_upload</v-icon>
        </v-btn>
        {{this.filetype}}
      </v-tooltip>
    </v-text-field>
        <div v-for="file in files" :key="file.id">
          <v-text-field 
            readonly 
            loading
            @click="$refs.upload.$el.lastChild.click()"
            @click:clear="remove(file)"
            v-model="file.name" 
            append-outer-icon="clear"
            @click:append-outer="remove(file)"
            class='ma-0'>
          <v-progress-linear
                  color="accent"
                  height="3"
                  slot="progress"
                  :value="file.progress"
            ></v-progress-linear>
            </v-text-field>
        </div>
      <v-btn small :disabled='getNumFilesWaitingForUpload() == 0' v-if="!$refs.upload || !$refs.upload.active" @click.prevent="$refs.upload.active = true">
        <v-icon small class="fa fa-arrow-up" aria-hidden="true"></v-icon>&nbsp;
        Start Upload
      </v-btn>
      <v-btn small v-else @click.prevent="stopUpload()">
        <v-icon small class="fa fa-stop" aria-hidden="true"></v-icon>&nbsp;
        Stop Upload
      </v-btn>
      <file-upload
      class="btn btn-primary"
      :custom-action="uploadFile"
      :extensions="filesuffix"
      :accept='filesuffix'
      :multiple="true"
      :size="sizelimit*1024*1024"
      v-model="files"
      @input-filter="inputFilter"
      @input-file="inputFile"
      :data="{type: filetype}"
      ref="upload">
    </file-upload>
  </v-container>
</template>

<script>
import FileUpload from 'vue-upload-component'
import axios from 'axios'
import store from '@/store'

export default {
  components: {
    FileUpload,
  },
  
  props: {
    'filetype': String, 
    'hinttext': String,
    'filesuffix': String,
    'sizelimit': {
      type: Number,
      default: -1
    }
  },
  
  data() {
    return {
      files: [],
      progress: 0,
      baseurl: '/prosit/api/upload.xsjs',
      file: null,
      filename: this.filetype,
      requesttoken: axios.CancelToken.source()
    }
  },

  methods: {
    remove(file) {
      this.$refs.upload.remove(file)
      this.checkNumUploadFiles()
    },
    stopUpload() {
      this.requesttoken.cancel();
      this.resetUpload()
    },
    resetUpload() {
      this.$refs.upload.active = false
    },
    checkNumUploadFiles() {
      var count = 0
      this.files.forEach(function (f) {
        if (f.success) {
          count += 1
        }
      })
      if (count == 0) {
        store.commit('reportUpload', {fileType: this.filetype, isSuccessful: false});
      }
    },
    getNumFilesWaitingForUpload() {
      var count = 0
      this.files.forEach(function (f) {
        if (!f.success) {
          count += 1
        }
      })
      return count
    },
    inputFilter(newFile, oldFile, prevent) {
      if (newFile && !oldFile) {
        // Before adding a file

        // Filter system files or hide files
        if (/(\/|^)(Thumbs\.db|desktop\.ini|\..+)$/.test(newFile.name)) {
          return prevent()
        }

        // Filter php html js file
        if (/\.(php5?|html?|jsx?)$/i.test(newFile.name)) {
          return prevent()
        }
      }
    },
    inputFile(newFile, oldFile) {
      if (newFile && !oldFile) {
        // add
        //console.log('add', newFile)
        store.commit('reportUpload', {fileType: this.filetype, isSuccessful: false});
      }
      
      if (newFile && oldFile) {
        // update
        //console.log('update', newFile, oldFile)
        
        // on upload abort
        if (!newFile.active && oldFile.active && newFile.progress != 100) {
          newFile.progress = 0
          newFile.success = false
        }
      }

      if (!newFile && oldFile) {
        // remove
        //console.log('remove', oldFile)
      }
    },
    async uploadFile(file) {
      this.requesttoken.cancel();
      file.progress = 0;

      this.file = file;
      this.filename = this.file.name;
      let formData = new FormData();
      formData.set('type', this.filetype);
      formData.append('file', this.file.file);
      this.requesttoken = axios.CancelToken.source();

      await axios.post(
        this.baseurl,
        formData,
        {
          params: { datasetId: this.$store.state.taskID },
          headers: { 'Content-Type': 'multipart/form-data' },
          onUploadProgress: ( progressEvent ) => {
            file.progress = parseInt( Math.round( ( progressEvent.loaded * 100 ) / progressEvent.total ) );
          },
          cancelToken: this.requesttoken.token
        }
      ).then(() => {
        // success
        if (this.getNumFilesWaitingForUpload() == 1) {
          store.commit('reportUpload', {fileType: this.filetype, isSuccessful: true});
        }
      }).catch(() =>{
        // error
        alert('Something went wrong. Please reload the page and try again (respectfully).')
        this.resetUpload()
      });
    },
  }
}
</script>
