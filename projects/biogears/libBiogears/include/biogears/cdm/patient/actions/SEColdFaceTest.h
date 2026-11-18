/**************************************************************************************
Copyright 2026 INESC TEC / FEUP-MEF (extension to BioGears v5 — Lucas Teixeira).
Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the License
at:
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under
the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.
**************************************************************************************/
//
//  SEColdFaceTest.h
//  ----------------
//  Action que dispara o reflexo trigemino-cardiaco (Cold Face Test) no
//  controlador nervoso. Modela:
//      - Aumento aditivo do sinal vagal eferente ao no SA (bradicardia)
//      - Aumento aditivo do sinal simpatico periferico (vasoconstricao)
//      - Sem modulacao do sinal simpatico-sinoatrial (preserva o caracter
//        non-baroreflex demonstrado por Khurana & Wu 2006).
//
//  Calibrado contra Heath & Downey 1990 (Clin Sci 78:139-147):
//      - Cold pack a 0-4 °C aplicado a regiao V1 do trigemeo
//      - Bradicardia esperada: -10 a -25% do HR basal
//      - Latencia ao pico: 30-60 s
//      - Vasoconstricao periferica com ΔSBP +5 a +15 mmHg
//
//  Parametros:
//      - Severity [0,1]: amplitude da modulacao (1.0 = calibracao Heath&Downey)
//      - Duration [s]:   duracao total do estimulo (default 90 s)
//                        Apos Duration, a fase temporal phi cai a zero e a
//                        modulacao torna-se inerte; a action permanece na
//                        collection ate ser explicitamente removida via
//                        SetSeverity(0) ou RemoveColdFaceTest().
//
//  Uso programatico (exemplo):
//      SEColdFaceTest cft;
//      cft.GetSeverity().SetValue(1.0);
//      cft.GetDuration().SetValue(90.0, TimeUnit::s);
//      bg->ProcessAction(cft);
//
//  NOTA: esta action NAO suporta serializacao XML (Marshall/UnMarshall) por
//        nao haver schema CDM associado. Funcionalidade restrita a uso
//        programatico em codigo C++. Para suporte XML, sera necessario
//        adicionar ColdFaceTestData ao XSD e regenerar Properties.hxx.
//=====================================================================*/

#pragma once
#include <biogears/cdm/patient/actions/SEPatientAction.h>
#include <biogears/cdm/enums/SEPatientActionsEnums.h>
#include <biogears/cdm/properties/SEScalarTime.h>

namespace biogears {
class SEScalar0To1;


class BIOGEARS_API SEColdFaceTest : public SEPatientAction {
public:
  SEColdFaceTest();
  virtual ~SEColdFaceTest() override;

  static constexpr const char* TypeTag() { return "SEColdFaceTest"; };
  const char* classname() const override { return TypeTag(); }

  virtual void Invalidate() override;

  virtual bool IsValid() const override;
  virtual bool IsActive() const override;

  virtual bool HasSeverity() const;
  virtual SEScalar0To1& GetSeverity();

  virtual bool HasDuration() const;
  virtual SEScalarTime& GetDuration();

  virtual void ToString(std::ostream& str) const override;

  bool operator==(const SEColdFaceTest& rhs) const;
  bool operator!=(const SEColdFaceTest& rhs) const;

protected:
  SEScalar0To1* m_Severity;
  SEScalarTime* m_Duration;
};
} // namespace biogears
