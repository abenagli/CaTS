//
// class ColorReader
//
// Custom reader for handling "color" tags extensions in GDML.
// -------------------------------------------------------------------------

#ifndef ColorReader_H
#define ColorReader_H 1

#include <map>
#include "G4GDMLReadStructure.hh"

class G4VisAttributes;

/// GDML reader for the color attributes

class ColorReader : public G4GDMLReadStructure
{

 public:

   ColorReader();
  ~ColorReader();

   void ExtensionRead(const xercesc::DOMElement* const element);
   void ColorRead(const xercesc::DOMElement* const element);

   G4VisAttributes* GetVisAttribute(const G4String& ref);

 protected:

   virtual void VolumeRead(const xercesc::DOMElement* const);

 private:

   std::map<G4String, G4VisAttributes*> fAttribs;
};

#endif
