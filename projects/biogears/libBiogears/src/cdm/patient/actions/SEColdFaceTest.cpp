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
#include <biogears/cdm/patient/actions/SEColdFaceTest.h>

#include <biogears/cdm/properties/SEScalar0To1.h>
#include <biogears/cdm/properties/SEScalarTime.h>

namespace biogears {

SEColdFaceTest::SEColdFaceTest()
  : SEPatientAction()
{
  m_Severity = nullptr;
  m_Duration = nullptr;
}
//-------------------------------------------------------------------------------
SEColdFaceTest::~SEColdFaceTest()
{
  Invalidate();
}
//-------------------------------------------------------------------------------
void SEColdFaceTest::Invalidate()
{
  SEPatientAction::Invalidate();
  SAFE_DELETE(m_Severity);
  SAFE_DELETE(m_Duration);
}
//-------------------------------------------------------------------------------
bool SEColdFaceTest::IsValid() const
{
  // Severity e mandatorio. Duration tem default na leitura do controlador (90 s).
  return SEPatientAction::IsValid() && HasSeverity();
}
//-------------------------------------------------------------------------------
bool SEColdFaceTest::IsActive() const
{
  return IsValid() ? !m_Severity->IsZero() : false;
}
//-------------------------------------------------------------------------------
bool SEColdFaceTest::HasSeverity() const
{
  return m_Severity == nullptr ? false : m_Severity->IsValid();
}
//-------------------------------------------------------------------------------
SEScalar0To1& SEColdFaceTest::GetSeverity()
{
  if (m_Severity == nullptr)
    m_Severity = new SEScalar0To1();
  return *m_Severity;
}
//-------------------------------------------------------------------------------
bool SEColdFaceTest::HasDuration() const
{
  return m_Duration == nullptr ? false : m_Duration->IsValid();
}
//-------------------------------------------------------------------------------
SEScalarTime& SEColdFaceTest::GetDuration()
{
  if (m_Duration == nullptr)
    m_Duration = new SEScalarTime();
  return *m_Duration;
}
//-------------------------------------------------------------------------------
void SEColdFaceTest::ToString(std::ostream& str) const
{
  str << "Patient Action : Cold Face Test (Trigeminocardiac Reflex)";
  if (HasComment())
    str << "\n\tComment: " << m_Comment;
  str << "\n\tSeverity: ";
  HasSeverity() ? str << *m_Severity : str << "Not Set";
  str << "\n\tDuration: ";
  HasDuration() ? str << *m_Duration : str << "Not Set (default 90 s)";
  str << std::flush;
}
//-------------------------------------------------------------------------------
bool SEColdFaceTest::operator==(const SEColdFaceTest& rhs) const
{
  bool equivilant = m_Comment == rhs.m_Comment;
  equivilant &= (m_Severity && rhs.m_Severity)
                  ? m_Severity->operator==(*rhs.m_Severity)
                  : m_Severity == rhs.m_Severity;
  equivilant &= (m_Duration && rhs.m_Duration)
                  ? m_Duration->operator==(*rhs.m_Duration)
                  : m_Duration == rhs.m_Duration;
  return equivilant;
}
//-------------------------------------------------------------------------------
bool SEColdFaceTest::operator!=(const SEColdFaceTest& rhs) const
{
  return !(*this == rhs);
}
} // namespace biogears
