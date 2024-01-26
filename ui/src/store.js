import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
  state: {
    taskTypes: ['CollisionEnergyAlignment', 'SpectralLibraryGeneration','MaxQuantRescoring'],//,'predictProperties',
    task: {
      type: "MaxQuantRescoring",
      hashId: false,
      taskId: false,
      tag: 'tmt',
      allFeatures: false,
      uploads: { 
        'msms.txt': false, 
        'RAW': false, 
        'FASTA': false, 
        'peptides.csv': false
      },
      options: { 
        'collisionEnergy': false,
        'protease': false,
        'missedCleavages': false,
        'oxidizedMethionine': false 
      },
      outputFormat: false,
      selectedIntensityModel: false,
      selectedIRTModel: false,
      searchEngine: false,
    },
  },
  mutations: {
    reportUpload (state, upload){
      state.task.uploads[upload.fileType] = upload.isSuccessful;
    },
    changeTask(state, taskType){
      state.task.type = taskType;
    },
    changeFormat(state, taskOutputFormat){
      state.task.outputFormat = taskOutputFormat;
    },
    changeIRTModel(state, irtModelName){
      state.task.selectedIRTModel = irtModelName;
    },
    changeIntensityModel(state, intensityModelName){
      state.task.selectedIntensityModel = intensityModelName;
    },
    changeTag(state, tag){
      state.task.tag = tag;
    },
    changeSearchEngine(state, searchEngine){
      state.task.searchEngine = searchEngine;
    },
    sethashId(state, hashId){
      state.task.hashId = hashId;
    },
    setTaskId(state, taskId){
      state.task.taskId = taskId;
    },
    changeOptions(state, options){
      for (const option in options) {
        if (state.task.options.hasOwnProperty(option)) {
          state.task.options[option] = options[option]
        }
      }
    },
    reset (state){
      state.task.type = "MaxQuantRescoring";
      for (const fileType in state.task.uploads) {
        if (state.task.uploads.hasOwnProperty(fileType)) {
          state.task.uploads[fileType] = false
        }
      }
      for (const fileType in state.task.options) {
        if (state.task.options.hasOwnProperty(fileType)) {
          state.task.options[fileType] = false
        }
      }
    }
  },
  getters: {
    submitObject: state => {
      let task = JSON.parse(JSON.stringify(state.task))
      return task
    },
    getTask: state => {
      let task = state.task
      return task
    },
  }
})
