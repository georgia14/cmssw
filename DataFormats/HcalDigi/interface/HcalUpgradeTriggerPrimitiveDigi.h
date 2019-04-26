#ifndef DIGIHCAL_HCALUPGRADETRIGGERPRIMITIVEDIGI_H
#define DIGIHCAL_HCALUPGRADETRIGGERPRIMITIVEDIGI_H

#include <ostream>
#include <vector>
#include "DataFormats/HcalDetId/interface/HcalTrigTowerDetId.h"
#include "DataFormats/HcalDigi/interface/HcalUpgradeTriggerPrimitiveSample.h"

class HcalUpgradeTriggerPrimitiveDigi {

 public:
  
  //------------------------------------------------------
  // Needed for consistency with SortedCollection
  //------------------------------------------------------
  
  typedef HcalTrigTowerDetId key_type;
  explicit HcalUpgradeTriggerPrimitiveDigi(const HcalTrigTowerDetId& id);
  const HcalTrigTowerDetId& id() const { return m_id; }

  //------------------------------------------------------
  // Constructor/Destructor
  //------------------------------------------------------
  
  HcalUpgradeTriggerPrimitiveDigi(); 
  ~HcalUpgradeTriggerPrimitiveDigi();

  //------------------------------------------------------
  // Set information.  MAXSAMPLES sets size limit.
  //------------------------------------------------------

  void setDepthData ( std::vector<int> d ) { m_depth_data = d; };

  void setTimingData (
		      std::vector<double> rise_avg,
		      std::vector<double> rise_rms,
		      std::vector<double> fall_avg,
		      std::vector<double> fall_rms
		      ) {
    m_rising_avg = rise_avg;
    m_rising_rms = rise_rms;
    m_falling_avg = fall_avg;
    m_falling_rms = fall_rms;
  };
  void setOOTData ( std::vector<int> d ) { m_oot_data = d; };

  void setSize      ( int  size );
  void setPresamples( int  presamples );
  void setZSInfo    ( bool unsuppressed, bool markAndPass);
  void setSample    ( int i, const HcalUpgradeTriggerPrimitiveSample& sample ) { m_data[i] = sample; }

  const std::vector<int>& getDepthData() const { return m_depth_data; };

  const std::vector<double>& getRisingAvg() const { return m_rising_avg; };
  const std::vector<double>& getRisingRMS() const { return m_rising_rms; };
  const std::vector<double>& getFallingAvg() const { return m_falling_avg; };
  const std::vector<double>& getFallingRMS() const { return m_falling_rms; };

  static const int MAXSAMPLES = 10;

  //------------------------------------------------------
  // Get the number of samples / presamples
  //------------------------------------------------------
  
  int size      () const { return (m_size           & 0xF); } // 
  int presamples() const { return (m_hcalPresamples & 0xF); } 
  
  //------------------------------------------------------
  // Get nformation about the ZS
  //------------------------------------------------------

  bool zsMarkAndPass () const { return m_hcalPresamples & 0x10; }
  bool zsUnsuppressed() const { return m_hcalPresamples & 0x20; }

  //------------------------------------------------------
  // Get information about individual samples
  //------------------------------------------------------
  
  // Access all stored samples

  const HcalUpgradeTriggerPrimitiveSample& operator[](int i) const { return m_data[i]; }
  const HcalUpgradeTriggerPrimitiveSample& sample    (int i) const { return m_data[i]; }
  
  // Access "sample of interest" directly
  
  const HcalUpgradeTriggerPrimitiveSample& t0() const { return m_data[presamples()]; }
  int SOI_fineGrain      () const { return t0().fineGrain      (); }
  int SOI_compressedEt   () const { return t0().compressedEt   (); }

  int SOI_depth_linear(int i) const { return m_depth_data[i - 1]; }

  double SOI_rising_avg(int i) const { return m_rising_avg[i - 1]; }
  double SOI_rising_rms(int i) const { return m_rising_rms[i - 1]; }
  double SOI_falling_avg(int i) const { return m_falling_avg[i - 1]; }
  double SOI_falling_rms(int i) const { return m_falling_rms[i - 1]; }
  int SOI_oot_linear(int i) const { return m_oot_data[i - 1]; }


 private:
  
  HcalTrigTowerDetId m_id;
  int m_size;
  int m_hcalPresamples;
  HcalUpgradeTriggerPrimitiveSample m_data [MAXSAMPLES];

  std::vector<int> m_depth_data;

  std::vector<int> m_oot_data;
  std::vector<double> m_rising_avg;
  std::vector<double> m_rising_rms;
  std::vector<double> m_falling_avg;
  std::vector<double> m_falling_rms;

};

#endif
