/**************************************************************************************
Copyright 2015 Applied Research Associates, Inc.
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the License
at:
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
**************************************************************************************/

// Include the various types you will be using in your code
#include <biogears/cdm/compartment/SECompartmentManager.h>
#include <biogears/cdm/engine/PhysiologyEngineTrack.h>
#include <biogears/cdm/patient/actions/SEAsthmaAttack.h>
#include <biogears/cdm/patient/actions/SESubstanceBolus.h>
#include <biogears/cdm/patient/actions/SESubstanceOralDose.h>
#include <biogears/cdm/patient/actions/SEExercise.h>
#include <biogears/cdm/properties/SEProperties.h>
#include <biogears/cdm/system/physiology/SEBloodChemistrySystem.h>
#include <biogears/cdm/system/physiology/SEEndocrineSystem.h>
#include <biogears/cdm/system/physiology/SEHepaticSystem.h>
#include <biogears/cdm/system/physiology/SEEnergySystem.h>
#include <biogears/cdm/system/physiology/SECardiovascularSystem.h>
#include <biogears/cdm/system/physiology/SERespiratorySystem.h>
#include <biogears/cdm/system/physiology/SETissueSystem.h>
#include <biogears/cdm/system/equipment/ElectroCardioGram/SEElectroCardioGram.h>
#include <biogears/cdm/engine/PhysiologyEngineConfiguration.h>
#include <biogears/cdm/system/physiology/SEDrugSystem.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <biogears/string/manipulation.h>
#include <biogears/engine/BioGearsPhysiologyEngine.h>
#include <fstream>

#include <sstream>

using namespace biogears;
//--------------------------------------------------------------------------------------------------
/// \brief
///
/// \details
//--------------------------------------------------------------------------------------------------
int HowToTesteP01()
{
  // Create the engine and load the patient
  Logger logger("HowToP01.log");
  std::unique_ptr<PhysiologyEngine> bg = CreateBioGearsEngine(&logger);
  

 
  auto& drm = bg->GetEngineTrack()->GetDataRequestManager();

  drm.SetSamplesPerSecond(250);
  
  drm.CreatePhysiologyDataRequest().Set("HeartRate", FrequencyUnit::Per_min);
  drm.CreatePhysiologyDataRequest().Set("SkinTemperature", TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("CoreTemperature", TemperatureUnit::C);
  drm.CreatePhysiologyDataRequest().Set("InspiratoryExpiratoryRatio");
  drm.CreatePhysiologyDataRequest().Set("CardiacOutput", VolumePerTimeUnit::mL_Per_min);
  drm.CreatePhysiologyDataRequest().Set("MeanArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("SystolicArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("DiastolicArterialPressure", PressureUnit::mmHg);
  drm.CreatePhysiologyDataRequest().Set("HemoglobinContent", MassUnit::g);
  drm.CreatePhysiologyDataRequest().Set("GlucagonSynthesisRate", AmountPerTimeUnit::umol_Per_s);
  drm.CreatePhysiologyDataRequest().Set("InsulinSynthesisRate", AmountPerTimeUnit::umol_Per_s);
  drm.CreatePhysiologyDataRequest().Set("LactateProductionRate", AmountPerTimeUnit::umol_Per_s);
  drm.CreatePhysiologyDataRequest().Set("FatigueLevel");
  drm.CreatePhysiologyDataRequest().Set("CreatinineProductionRate", AmountPerTimeUnit::umol_Per_s);
  drm.CreatePhysiologyDataRequest().Set("ChlorideLostToSweat", MassUnit::g);
  drm.CreatePhysiologyDataRequest().Set("SodiumLostToSweat", MassUnit::g);
  drm.CreatePhysiologyDataRequest().Set("HepaticGluconeogenesisRate", MassPerTimeUnit::mg_Per_min);
  drm.CreatePhysiologyDataRequest().Set("OxygenConsumptionRate", VolumePerTimeUnit::mL_Per_s);

  drm.SetResultsFilename("HowToP01.csv");   // ← Ficheiro fisiológico principal

  if (!bg->InitializeEngine("patients/Teste01.xml")) {
    bg->GetLogger()->Error("Could not load state, check the error");
    return 1;
  }

  double dt = 1.0 / 250.0;

   // --- Cria ficheiro separado para ECG ---
  std::ofstream ecgCSV("RunningP01.csv");
  ecgCSV << "Time_s,"
         << "Lead3_mV\n";

  // --- Simula baseline (2 minutos = 600s) ---
  bg->GetLogger()->Info("Simulando baseline (2min)...");
  for (int t = 0; t < 120/dt; ++t) {
    bg->AdvanceModelTime(dt, TimeUnit::s);

    const auto& ecg = *bg->GetElectroCardioGram();
    double time_s = bg->GetSimulationTime(TimeUnit::s);

    ecgCSV << time_s << ","
           << ecg.GetLead3ElectricPotential(ElectricPotentialUnit::mV) << "\n";
  }


  // Running starts

  SEExercise exercise;
  SEExercise::SERunning running;
  
  running.SpeedRun.SetValue(2.4, LengthPerTimeUnit::m_Per_s); // define velocidade de corrida (exemplo)
  running.InclineRun.SetValue(0.0); // corrida em plano
  exercise.SetRunningExercise(running);
  bg->GetLogger()->Info("The patient starts running ");
  bg->ProcessAction(exercise); // garante que HasExercise() continua true

  for (int t = 0; t < 180/dt; ++t) {
    bg->AdvanceModelTime(dt, TimeUnit::s);
    const auto& ecg = *bg->GetElectroCardioGram();
    double time_s = bg->GetSimulationTime(TimeUnit::s);

    ecgCSV << time_s << ","
           << ecg.GetLead3ElectricPotential(ElectricPotentialUnit::mV) << "\n";
  }

  // REMOVER o exercício para entrar em recuperação
  bg->GetLogger()->Info("The patient stops running and enters recovery");

  SEExercise recovery; // Tipo NONE por padrão
  bg->ProcessAction(recovery); // Processa ação sem exercício

  for (int t = 0; t < 300/dt; ++t) {
    
      bg->AdvanceModelTime(dt, TimeUnit::s);
      const auto& ecg = *bg->GetElectroCardioGram();
      double time_s = bg->GetSimulationTime(TimeUnit::s);
      
      ecgCSV << time_s << ","
            << ecg.GetLead3ElectricPotential(ElectricPotentialUnit::mV) << "\n";
  }

  bg->GetLogger()->Info("Finished");
  return 0;
  
}

int main ( int argc, char* argv[] ) {
  return HowToTesteP01();
}
