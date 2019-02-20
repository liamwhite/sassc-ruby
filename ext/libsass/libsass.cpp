#include "libsass.h"

namespace Sass {

  /* the conversion matrix can be readed the following way */
  /* if you go down, the factor is for the numerator (multiply) */
  /* if you go right, the factor is for the denominator (divide) */
  /* and yes, we actually use both, not sure why, but why not!? */

  const double size_conversion_factors[6][6] =
  {
             /*  in         cm         pc         mm         pt         px        */
    /* in   */ { 1,         2.54,      6,         25.4,      72,        96,       },
    /* cm   */ { 1.0/2.54,  1,         6.0/2.54,  10,        72.0/2.54, 96.0/2.54 },
    /* pc   */ { 1.0/6.0,   2.54/6.0,  1,         25.4/6.0,  72.0/6.0,  96.0/6.0  },
    /* mm   */ { 1.0/25.4,  1.0/10.0,  6.0/25.4,  1,         72.0/25.4, 96.0/25.4 },
    /* pt   */ { 1.0/72.0,  2.54/72.0, 6.0/72.0,  25.4/72.0, 1,         96.0/72.0 },
    /* px   */ { 1.0/96.0,  2.54/96.0, 6.0/96.0,  25.4/96.0, 72.0/96.0, 1,        }
  };

  const double angle_conversion_factors[4][4] =
  {
             /*  deg        grad       rad        turn      */
    /* deg  */ { 1,         40.0/36.0, PI/180.0,  1.0/360.0 },
    /* grad */ { 36.0/40.0, 1,         PI/200.0,  1.0/400.0 },
    /* rad  */ { 180.0/PI,  200.0/PI,  1,         0.5/PI    },
    /* turn */ { 360.0,     400.0,     2.0*PI,    1         }
  };

  const double time_conversion_factors[2][2] =
  {
             /*  s          ms        */
    /* s    */ { 1,         1000.0    },
    /* ms   */ { 1/1000.0,  1         }
  };
  const double frequency_conversion_factors[2][2] =
  {
             /*  Hz         kHz       */
    /* Hz   */ { 1,         1/1000.0  },
    /* kHz  */ { 1000.0,    1         }
  };
  const double resolution_conversion_factors[3][3] =
  {
             /*  dpi        dpcm       dppx     */
    /* dpi  */ { 1,         1/2.54,    1/96.0   },
    /* dpcm */ { 2.54,      1,         2.54/96  },
    /* dppx */ { 96,        96/2.54,   1        }
  };

  UnitClass get_unit_type(UnitType unit)
  {
    switch (unit & 0xFF00)
    {
      case UnitClass::LENGTH:      return UnitClass::LENGTH;
      case UnitClass::ANGLE:       return UnitClass::ANGLE;
      case UnitClass::TIME:        return UnitClass::TIME;
      case UnitClass::FREQUENCY:   return UnitClass::FREQUENCY;
      case UnitClass::RESOLUTION:  return UnitClass::RESOLUTION;
      default:                     return UnitClass::INCOMMENSURABLE;
    }
  };

  std::string get_unit_class(UnitType unit)
  {
    switch (unit & 0xFF00)
    {
      case UnitClass::LENGTH:      return "LENGTH";
      case UnitClass::ANGLE:       return "ANGLE";
      case UnitClass::TIME:        return "TIME";
      case UnitClass::FREQUENCY:   return "FREQUENCY";
      case UnitClass::RESOLUTION:  return "RESOLUTION";
      default:                     return "INCOMMENSURABLE";
    }
  };

  UnitType get_main_unit(const UnitClass unit)
  {
    switch (unit)
    {
      case UnitClass::LENGTH:      return UnitType::PX;
      case UnitClass::ANGLE:       return UnitType::DEG;
      case UnitClass::TIME:        return UnitType::SEC;
      case UnitClass::FREQUENCY:   return UnitType::HERTZ;
      case UnitClass::RESOLUTION:  return UnitType::DPI;
      default:                     return UnitType::UNKNOWN;
    }
  };

  UnitType string_to_unit(const std::string& s)
  {
    // size units
    if      (s == "px")   return UnitType::PX;
    else if (s == "pt")   return UnitType::PT;
    else if (s == "pc")   return UnitType::PC;
    else if (s == "mm")   return UnitType::MM;
    else if (s == "cm")   return UnitType::CM;
    else if (s == "in")   return UnitType::IN;
    // angle units
    else if (s == "deg")  return UnitType::DEG;
    else if (s == "grad") return UnitType::GRAD;
    else if (s == "rad")  return UnitType::RAD;
    else if (s == "turn") return UnitType::TURN;
    // time units
    else if (s == "s")    return UnitType::SEC;
    else if (s == "ms")   return UnitType::MSEC;
    // frequency units
    else if (s == "Hz")   return UnitType::HERTZ;
    else if (s == "kHz")  return UnitType::KHERTZ;
    // resolutions units
    else if (s == "dpi")  return UnitType::DPI;
    else if (s == "dpcm") return UnitType::DPCM;
    else if (s == "dppx") return UnitType::DPPX;
    // for unknown units
    else return UnitType::UNKNOWN;
  }

  const char* unit_to_string(UnitType unit)
  {
    switch (unit) {
      // size units
      case UnitType::PX:      return "px";
      case UnitType::PT:      return "pt";
      case UnitType::PC:      return "pc";
      case UnitType::MM:      return "mm";
      case UnitType::CM:      return "cm";
      case UnitType::IN:      return "in";
      // angle units
      case UnitType::DEG:     return "deg";
      case UnitType::GRAD:    return "grad";
      case UnitType::RAD:     return "rad";
      case UnitType::TURN:    return "turn";
      // time units
      case UnitType::SEC:     return "s";
      case UnitType::MSEC:    return "ms";
      // frequency units
      case UnitType::HERTZ:   return "Hz";
      case UnitType::KHERTZ:  return "kHz";
      // resolutions units
      case UnitType::DPI:     return "dpi";
      case UnitType::DPCM:    return "dpcm";
      case UnitType::DPPX:    return "dppx";
      // for unknown units
      default:                return "";
    }
  }

  std::string unit_to_class(const std::string& s)
  {
    if      (s == "px")   return "LENGTH";
    else if (s == "pt")   return "LENGTH";
    else if (s == "pc")   return "LENGTH";
    else if (s == "mm")   return "LENGTH";
    else if (s == "cm")   return "LENGTH";
    else if (s == "in")   return "LENGTH";
    // angle units
    else if (s == "deg")  return "ANGLE";
    else if (s == "grad") return "ANGLE";
    else if (s == "rad")  return "ANGLE";
    else if (s == "turn") return "ANGLE";
    // time units
    else if (s == "s")    return "TIME";
    else if (s == "ms")   return "TIME";
    // frequency units
    else if (s == "Hz")   return "FREQUENCY";
    else if (s == "kHz")  return "FREQUENCY";
    // resolutions units
    else if (s == "dpi")  return "RESOLUTION";
    else if (s == "dpcm") return "RESOLUTION";
    else if (s == "dppx") return "RESOLUTION";
    // for unknown units
    return "CUSTOM:" + s;
  }

  // throws incompatibleUnits exceptions
  double conversion_factor(const std::string& s1, const std::string& s2)
  {
    // assert for same units
    if (s1 == s2) return 1;
    // get unit enum from string
    UnitType u1 = string_to_unit(s1);
    UnitType u2 = string_to_unit(s2);
    // query unit group types
    UnitClass t1 = get_unit_type(u1);
    UnitClass t2 = get_unit_type(u2);
    // return the conversion factor
    return conversion_factor(u1, u2, t1, t2);
  }

  // throws incompatibleUnits exceptions
  double conversion_factor(UnitType u1, UnitType u2, UnitClass t1, UnitClass t2)
  {
    // can't convert between groups
    if (t1 != t2) return 0;
    // get absolute offset
    // used for array acces
    size_t i1 = u1 - t1;
    size_t i2 = u2 - t2;
    // process known units
    switch (t1) {
      case LENGTH:
        return size_conversion_factors[i1][i2];
      case ANGLE:
        return angle_conversion_factors[i1][i2];
      case TIME:
        return time_conversion_factors[i1][i2];
      case FREQUENCY:
        return frequency_conversion_factors[i1][i2];
      case RESOLUTION:
        return resolution_conversion_factors[i1][i2];
      case INCOMMENSURABLE:
        return 0;
    }
    // fallback
    return 0;
  }

  double convert_units(const std::string& lhs, const std::string& rhs, int& lhsexp, int& rhsexp)
  {
    double f = 0;
    // do not convert same ones
    if (lhs == rhs) return 0;
    // skip already canceled out unit
    if (lhsexp == 0) return 0;
    if (rhsexp == 0) return 0;
    // check if it can be converted
    UnitType ulhs = string_to_unit(lhs);
    UnitType urhs = string_to_unit(rhs);
    // skip units we cannot convert
    if (ulhs == UNKNOWN) return 0;
    if (urhs == UNKNOWN) return 0;
    // query unit group types
    UnitClass clhs = get_unit_type(ulhs);
    UnitClass crhs = get_unit_type(urhs);
    // skip units we cannot convert
    if (clhs != crhs) return 0;
    // if right denominator is bigger than lhs, we want to keep it in rhs unit
    if (rhsexp < 0 && lhsexp > 0 && - rhsexp > lhsexp) {
      // get the conversion factor for units
      f = conversion_factor(urhs, ulhs, clhs, crhs);
      // left hand side has been consumned
      f = std::pow(f, lhsexp);
      rhsexp += lhsexp;
      lhsexp = 0;
    }
    else {
      // get the conversion factor for units
      f = conversion_factor(ulhs, urhs, clhs, crhs);
      // right hand side has been consumned
      f = std::pow(f, rhsexp);
      lhsexp += rhsexp;
      rhsexp = 0;
    }
    return f;
  }

  bool Units::operator< (const Units& rhs) const
  {
    return (numerators < rhs.numerators) &&
           (denominators < rhs.denominators);
  }
  bool Units::operator== (const Units& rhs) const
  {
    return (numerators == rhs.numerators) &&
           (denominators == rhs.denominators);
  }
  bool Units::operator!= (const Units& rhs) const
  {
    return ! (*this == rhs);
  }

  double Units::normalize()
  {

    size_t iL = numerators.size();
    size_t nL = denominators.size();

    // the final conversion factor
    double factor = 1;

    for (size_t i = 0; i < iL; i++) {
      std::string &lhs = numerators[i];
      UnitType ulhs = string_to_unit(lhs);
      if (ulhs == UNKNOWN) continue;
      UnitClass clhs = get_unit_type(ulhs);
      UnitType umain = get_main_unit(clhs);
      if (ulhs == umain) continue;
      double f(conversion_factor(umain, ulhs, clhs, clhs));
      if (f == 0) throw std::runtime_error("INVALID");
      numerators[i] = unit_to_string(umain);
      factor /= f;
    }

    for (size_t n = 0; n < nL; n++) {
      std::string &rhs = denominators[n];
      UnitType urhs = string_to_unit(rhs);
      if (urhs == UNKNOWN) continue;
      UnitClass crhs = get_unit_type(urhs);
      UnitType umain = get_main_unit(crhs);
      if (urhs == umain) continue;
      double f(conversion_factor(umain, urhs, crhs, crhs));
      if (f == 0) throw std::runtime_error("INVALID");
      denominators[n] = unit_to_string(umain);
      factor /= f;
    }

    std::sort (numerators.begin(), numerators.end());
    std::sort (denominators.begin(), denominators.end());

    // return for conversion
    return factor;
  }

  double Units::reduce()
  {

    size_t iL = numerators.size();
    size_t nL = denominators.size();

    // have less than two units?
    if (iL + nL < 2) return 1;

    // first make sure same units cancel each other out
    // it seems that a map table will fit nicely to do this
    // we basically construct exponents for each unit
    // has the advantage that they will be pre-sorted
    std::map<std::string, int> exponents;

    // initialize by summing up occurences in unit vectors
    // this will already cancel out equivalent units (e.q. px/px)
    for (size_t i = 0; i < iL; i ++) exponents[numerators[i]] += 1;
    for (size_t n = 0; n < nL; n ++) exponents[denominators[n]] -= 1;

    // the final conversion factor
    double factor = 1;

    // convert between compatible units
    for (size_t i = 0; i < iL; i++) {
      for (size_t n = 0; n < nL; n++) {
        std::string &lhs = numerators[i], &rhs = denominators[n];
        int &lhsexp = exponents[lhs], &rhsexp = exponents[rhs];
        double f(convert_units(lhs, rhs, lhsexp, rhsexp));
        if (f == 0) continue;
        factor /= f;
      }
    }

    // now we can build up the new unit arrays
    numerators.clear();
    denominators.clear();

    // recreate sorted units vectors
    for (auto exp : exponents) {
      int &exponent = exp.second;
      while (exponent > 0 && exponent --)
        numerators.push_back(exp.first);
      while (exponent < 0 && exponent ++)
        denominators.push_back(exp.first);
    }

    // return for conversion
    return factor;

  }

  std::string Units::unit() const
  {
    std::string u;
    size_t iL = numerators.size();
    size_t nL = denominators.size();
    for (size_t i = 0; i < iL; i += 1) {
      if (i) u += '*';
      u += numerators[i];
    }
    if (nL != 0) u += '/';
    for (size_t n = 0; n < nL; n += 1) {
      if (n) u += '*';
      u += denominators[n];
    }
    return u;
  }

  bool Units::is_unitless() const
  {
    return numerators.empty() &&
           denominators.empty();
  }

  bool Units::is_valid_css_unit() const
  {
    return numerators.size() <= 1 &&
           denominators.size() == 0;
  }

  // this does not cover all cases (multiple prefered units)
  double Units::convert_factor(const Units& r) const
  {

    std::vector<std::string> miss_nums(0);
    std::vector<std::string> miss_dens(0);
    // create copy since we need these for state keeping
    std::vector<std::string> r_nums(r.numerators);
    std::vector<std::string> r_dens(r.denominators);

    auto l_num_it = numerators.begin();
    auto l_num_end = numerators.end();

    bool l_unitless = is_unitless();
    auto r_unitless = r.is_unitless();

    // overall conversion
    double factor = 1;

    // process all left numerators
    while (l_num_it != l_num_end)
    {
      // get and increment afterwards
      const std::string l_num = *(l_num_it ++);

      auto r_num_it = r_nums.begin(), r_num_end = r_nums.end();

      bool found = false;
      // search for compatible numerator
      while (r_num_it != r_num_end)
      {
        // get and increment afterwards
        const std::string r_num = *(r_num_it);
        // get possible conversion factor for units
        double conversion = conversion_factor(l_num, r_num);
        // skip incompatible numerator
        if (conversion == 0) {
          ++ r_num_it;
          continue;
        }
        // apply to global factor
        factor *= conversion;
        // remove item from vector
        r_nums.erase(r_num_it);
        // found numerator
        found = true;
        break;
      }
      // maybe we did not find any
      // left numerator is leftover
      if (!found) miss_nums.push_back(l_num);
    }

    auto l_den_it = denominators.begin();
    auto l_den_end = denominators.end();

    // process all left denominators
    while (l_den_it != l_den_end)
    {
      // get and increment afterwards
      const std::string l_den = *(l_den_it ++);

      auto r_den_it = r_dens.begin();
      auto r_den_end = r_dens.end();

      bool found = false;
      // search for compatible denominator
      while (r_den_it != r_den_end)
      {
        // get and increment afterwards
        const std::string r_den = *(r_den_it);
        // get possible converstion factor for units
        double conversion = conversion_factor(l_den, r_den);
        // skip incompatible denominator
        if (conversion == 0) {
          ++ r_den_it;
          continue;
        }
        // apply to global factor
        factor /= conversion;
        // remove item from vector
        r_dens.erase(r_den_it);
        // found denominator
        found = true;
        break;
      }
      // maybe we did not find any
      // left denominator is leftover
      if (!found) miss_dens.push_back(l_den);
    }

    // check left-overs (ToDo: might cancel out?)
    if (miss_nums.size() > 0 && !r_unitless) {
      throw Exception::IncompatibleUnits(r, *this);
    }
    else if (miss_dens.size() > 0 && !r_unitless) {
      throw Exception::IncompatibleUnits(r, *this);
    }
    else if (r_nums.size() > 0 && !l_unitless) {
      throw Exception::IncompatibleUnits(r, *this);
    }
    else if (r_dens.size() > 0 && !l_unitless) {
      throw Exception::IncompatibleUnits(r, *this);
    }

    return factor;
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

    Remove_Placeholders::Remove_Placeholders()
    { }

    void Remove_Placeholders::operator()(Block_Ptr b) {
        for (size_t i = 0, L = b->length(); i < L; ++i) {
            Statement_Ptr st = b->at(i);
            st->perform(this);
        }
    }

    Selector_List_Ptr Remove_Placeholders::remove_placeholders(Selector_List_Ptr sl)
    {
      Selector_List_Ptr new_sl = SASS_MEMORY_NEW(Selector_List, sl->pstate());

      for (size_t i = 0, L = sl->length(); i < L; ++i) {
          if (!sl->at(i)->contains_placeholder()) {
              new_sl->append(sl->at(i));
          }
      }

      return new_sl;

    }


    void Remove_Placeholders::operator()(Ruleset_Ptr r) {
        // Create a new selector group without placeholders
        Selector_List_Obj sl = Cast<Selector_List>(r->selector());

        if (sl) {
          // Set the new placeholder selector list
          r->selector(remove_placeholders(sl));
          // Remove placeholders in wrapped selectors
          for (Complex_Selector_Obj cs : sl->elements()) {
            while (cs) {
              if (cs->head()) {
                for (Simple_Selector_Obj& ss : cs->head()->elements()) {
                  if (Wrapped_Selector_Ptr ws = Cast<Wrapped_Selector>(ss)) {
                    if (Selector_List_Ptr wsl = Cast<Selector_List>(ws->selector())) {
                      Selector_List_Ptr clean = remove_placeholders(wsl);
                      // also clean superflous parent selectors
                      // probably not really the correct place
                      clean->remove_parent_selectors();
                      ws->selector(clean);
                    }
                  }
                }
              }
              cs = cs->tail();
            }
          }
        }

        // Iterate into child blocks
        Block_Obj b = r->block();

        for (size_t i = 0, L = b->length(); i < L; ++i) {
            if (b->at(i)) {
                Statement_Obj st = b->at(i);
                st->perform(this);
            }
        }
    }

    void Remove_Placeholders::operator()(Media_Block_Ptr m) {
        operator()(m->block());
    }
    void Remove_Placeholders::operator()(Supports_Block_Ptr m) {
        operator()(m->block());
    }

    void Remove_Placeholders::operator()(Directive_Ptr a) {
        if (a->block()) a->block()->perform(this);
    }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  void Subset_Map::put(const Compound_Selector_Obj& sel, const SubSetMapPair& value)
  {
    if (sel->empty()) throw std::runtime_error("internal error: subset map keys may not be empty");
    size_t index = values_.size();
    values_.push_back(value);
    for (size_t i = 0, S = sel->length(); i < S; ++i)
    {
      hash_[(*sel)[i]].push_back(std::make_pair(sel, index));
    }
  }

  std::vector<SubSetMapPair> Subset_Map::get_kv(const Compound_Selector_Obj& sel)
  {
    SimpleSelectorDict dict(sel->begin(), sel->end()); // XXX Set
    std::vector<size_t> indices;
    for (size_t i = 0, S = sel->length(); i < S; ++i) {
      if (!hash_.count((*sel)[i])) {
        continue;
      }
      const std::vector<std::pair<Compound_Selector_Obj, size_t> >& subsets = hash_[(*sel)[i]];
      for (const std::pair<Compound_Selector_Obj, size_t>& item : subsets) {
        bool include = true;
        for (const Simple_Selector_Obj& it : item.first->elements()) {
          auto found = dict.find(it);
          if (found == dict.end()) {
            include = false;
            break;
          }
        }
        if (include) indices.push_back(item.second);
      }
    }
    sort(indices.begin(), indices.end());
    std::vector<size_t>::iterator indices_end = unique(indices.begin(), indices.end());
    indices.resize(distance(indices.begin(), indices_end));

    std::vector<SubSetMapPair> results;
    for (size_t i = 0, S = indices.size(); i < S; ++i) {
      results.push_back(values_[indices[i]]);
    }
    return results;
  }

  std::vector<SubSetMapPair> Subset_Map::get_v(const Compound_Selector_Obj& sel)
  {
    return get_kv(sel);
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.

#include <iomanip>


namespace Sass {

  Inspect::Inspect(const Emitter& emi)
  : Emitter(emi)
  { }
  Inspect::~Inspect() { }

  // statements
  void Inspect::operator()(Block_Ptr block)
  {
    if (!block->is_root()) {
      add_open_mapping(block);
      append_scope_opener();
    }
    if (output_style() == NESTED) indentation += block->tabs();
    for (size_t i = 0, L = block->length(); i < L; ++i) {
      (*block)[i]->perform(this);
    }
    if (output_style() == NESTED) indentation -= block->tabs();
    if (!block->is_root()) {
      append_scope_closer();
      add_close_mapping(block);
    }

  }

  void Inspect::operator()(Ruleset_Ptr ruleset)
  {
    if (ruleset->selector()) {
      ruleset->selector()->perform(this);
    }
    if (ruleset->block()) {
      ruleset->block()->perform(this);
    }
  }

  void Inspect::operator()(Keyframe_Rule_Ptr rule)
  {
    if (rule->name()) rule->name()->perform(this);
    if (rule->block()) rule->block()->perform(this);
  }

  void Inspect::operator()(Bubble_Ptr bubble)
  {
    append_indentation();
    append_token("::BUBBLE", bubble);
    append_scope_opener();
    bubble->node()->perform(this);
    append_scope_closer();
  }

  void Inspect::operator()(Media_Block_Ptr media_block)
  {
    append_indentation();
    append_token("@media", media_block);
    append_mandatory_space();
    in_media_block = true;
    media_block->media_queries()->perform(this);
    in_media_block = false;
    media_block->block()->perform(this);
  }

  void Inspect::operator()(Supports_Block_Ptr feature_block)
  {
    append_indentation();
    append_token("@supports", feature_block);
    append_mandatory_space();
    feature_block->condition()->perform(this);
    feature_block->block()->perform(this);
  }

  void Inspect::operator()(At_Root_Block_Ptr at_root_block)
  {
    append_indentation();
    append_token("@at-root ", at_root_block);
    append_mandatory_space();
    if(at_root_block->expression()) at_root_block->expression()->perform(this);
    if(at_root_block->block()) at_root_block->block()->perform(this);
  }

  void Inspect::operator()(Directive_Ptr at_rule)
  {
    append_indentation();
    append_token(at_rule->keyword(), at_rule);
    if (at_rule->selector()) {
      append_mandatory_space();
      bool was_wrapped = in_wrapped;
      in_wrapped = true;
      at_rule->selector()->perform(this);
      in_wrapped = was_wrapped;
    }
    if (at_rule->value()) {
      append_mandatory_space();
      at_rule->value()->perform(this);
    }
    if (at_rule->block()) {
      at_rule->block()->perform(this);
    }
    else {
      append_delimiter();
    }
  }

  void Inspect::operator()(Declaration_Ptr dec)
  {
    if (dec->value()->concrete_type() == Expression::NULL_VAL) return;
    bool was_decl = in_declaration;
    in_declaration = true;
    LOCAL_FLAG(in_custom_property, dec->is_custom_property());

    if (output_style() == NESTED)
      indentation += dec->tabs();
    append_indentation();
    if (dec->property())
      dec->property()->perform(this);
    append_colon_separator();

    if (dec->value()->concrete_type() == Expression::SELECTOR) {
      Listize listize;
      Expression_Obj ls = dec->value()->perform(&listize);
      ls->perform(this);
    } else {
      dec->value()->perform(this);
    }

    if (dec->is_important()) {
      append_optional_space();
      append_string("!important");
    }
    append_delimiter();
    if (output_style() == NESTED)
      indentation -= dec->tabs();
    in_declaration = was_decl;
  }

  void Inspect::operator()(Assignment_Ptr assn)
  {
    append_token(assn->variable(), assn);
    append_colon_separator();
    assn->value()->perform(this);
    if (assn->is_default()) {
      append_optional_space();
      append_string("!default");
    }
    append_delimiter();
  }

  void Inspect::operator()(Import_Ptr import)
  {
    if (!import->urls().empty()) {
      append_token("@import", import);
      append_mandatory_space();

      import->urls().front()->perform(this);
      if (import->urls().size() == 1) {
        if (import->import_queries()) {
          append_mandatory_space();
          import->import_queries()->perform(this);
        }
      }
      append_delimiter();
      for (size_t i = 1, S = import->urls().size(); i < S; ++i) {
        append_mandatory_linefeed();
        append_token("@import", import);
        append_mandatory_space();

        import->urls()[i]->perform(this);
        if (import->urls().size() - 1 == i) {
          if (import->import_queries()) {
            append_mandatory_space();
            import->import_queries()->perform(this);
          }
        }
        append_delimiter();
      }
    }
  }

  void Inspect::operator()(Import_Stub_Ptr import)
  {
    append_indentation();
    append_token("@import", import);
    append_mandatory_space();
    append_string(import->imp_path());
    append_delimiter();
  }

  void Inspect::operator()(Warning_Ptr warning)
  {
    append_indentation();
    append_token("@warn", warning);
    append_mandatory_space();
    warning->message()->perform(this);
    append_delimiter();
  }

  void Inspect::operator()(Error_Ptr error)
  {
    append_indentation();
    append_token("@error", error);
    append_mandatory_space();
    error->message()->perform(this);
    append_delimiter();
  }

  void Inspect::operator()(Debug_Ptr debug)
  {
    append_indentation();
    append_token("@debug", debug);
    append_mandatory_space();
    debug->value()->perform(this);
    append_delimiter();
  }

  void Inspect::operator()(Comment_Ptr comment)
  {
    in_comment = true;
    comment->text()->perform(this);
    in_comment = false;
  }

  void Inspect::operator()(If_Ptr cond)
  {
    append_indentation();
    append_token("@if", cond);
    append_mandatory_space();
    cond->predicate()->perform(this);
    cond->block()->perform(this);
    if (cond->alternative()) {
      append_optional_linefeed();
      append_indentation();
      append_string("else");
      cond->alternative()->perform(this);
    }
  }

  void Inspect::operator()(For_Ptr loop)
  {
    append_indentation();
    append_token("@for", loop);
    append_mandatory_space();
    append_string(loop->variable());
    append_string(" from ");
    loop->lower_bound()->perform(this);
    append_string(loop->is_inclusive() ? " through " : " to ");
    loop->upper_bound()->perform(this);
    loop->block()->perform(this);
  }

  void Inspect::operator()(Each_Ptr loop)
  {
    append_indentation();
    append_token("@each", loop);
    append_mandatory_space();
    append_string(loop->variables()[0]);
    for (size_t i = 1, L = loop->variables().size(); i < L; ++i) {
      append_comma_separator();
      append_string(loop->variables()[i]);
    }
    append_string(" in ");
    loop->list()->perform(this);
    loop->block()->perform(this);
  }

  void Inspect::operator()(While_Ptr loop)
  {
    append_indentation();
    append_token("@while", loop);
    append_mandatory_space();
    loop->predicate()->perform(this);
    loop->block()->perform(this);
  }

  void Inspect::operator()(Return_Ptr ret)
  {
    append_indentation();
    append_token("@return", ret);
    append_mandatory_space();
    ret->value()->perform(this);
    append_delimiter();
  }

  void Inspect::operator()(Extension_Ptr extend)
  {
    append_indentation();
    append_token("@extend", extend);
    append_mandatory_space();
    extend->selector()->perform(this);
    append_delimiter();
  }

  void Inspect::operator()(Definition_Ptr def)
  {
    append_indentation();
    if (def->type() == Definition::MIXIN) {
      append_token("@mixin", def);
      append_mandatory_space();
    } else {
      append_token("@function", def);
      append_mandatory_space();
    }
    append_string(def->name());
    def->parameters()->perform(this);
    def->block()->perform(this);
  }

  void Inspect::operator()(Mixin_Call_Ptr call)
  {
    append_indentation();
    append_token("@include", call);
    append_mandatory_space();
    append_string(call->name());
    if (call->arguments()) {
      call->arguments()->perform(this);
    }
    if (call->block()) {
      append_optional_space();
      call->block()->perform(this);
    }
    if (!call->block()) append_delimiter();
  }

  void Inspect::operator()(Content_Ptr content)
  {
    append_indentation();
    append_token("@content", content);
    append_delimiter();
  }

  void Inspect::operator()(Map_Ptr map)
  {
    if (output_style() == TO_SASS && map->empty()) {
      append_string("()");
      return;
    }
    if (map->empty()) return;
    if (map->is_invisible()) return;
    bool items_output = false;
    append_string("(");
    for (auto key : map->keys()) {
      if (items_output) append_comma_separator();
      key->perform(this);
      append_colon_separator();
      LOCAL_FLAG(in_space_array, true);
      LOCAL_FLAG(in_comma_array, true);
      map->at(key)->perform(this);
      items_output = true;
    }
    append_string(")");
  }

  std::string Inspect::lbracket(List_Ptr list) {
    return list->is_bracketed() ? "[" : "(";
  }

  std::string Inspect::rbracket(List_Ptr list) {
    return list->is_bracketed() ? "]" : ")";
  }

  void Inspect::operator()(List_Ptr list)
  {
    if (list->empty() && (output_style() == TO_SASS || list->is_bracketed())) {
      append_string(lbracket(list));
      append_string(rbracket(list));
      return;
    }
    std::string sep(list->separator() == SASS_SPACE ? " " : ",");
    if ((output_style() != COMPRESSED) && sep == ",") sep += " ";
    else if (in_media_block && sep != " ") sep += " "; // verified
    if (list->empty()) return;
    bool items_output = false;

    bool was_space_array = in_space_array;
    bool was_comma_array = in_comma_array;
    // if the list is bracketed, always include the left bracket
    if (list->is_bracketed()) {
      append_string(lbracket(list));
    }
    // probably ruby sass eqivalent of element_needs_parens
    else if (output_style() == TO_SASS &&
        list->length() == 1 &&
        !list->from_selector() &&
        !Cast<List>(list->at(0)) &&
        !Cast<Selector_List>(list->at(0))
    ) {
      append_string(lbracket(list));
    }
    else if (!in_declaration && (list->separator() == SASS_HASH ||
        (list->separator() == SASS_SPACE && in_space_array) ||
        (list->separator() == SASS_COMMA && in_comma_array)
    )) {
      append_string(lbracket(list));
    }

    if (list->separator() == SASS_SPACE) in_space_array = true;
    else if (list->separator() == SASS_COMMA) in_comma_array = true;

    for (size_t i = 0, L = list->size(); i < L; ++i) {
      if (list->separator() == SASS_HASH)
      { sep[0] = i % 2 ? ':' : ','; }
      Expression_Obj list_item = list->at(i);
      if (output_style() != TO_SASS) {
        if (list_item->is_invisible()) {
          // this fixes an issue with "" in a list
          if (!Cast<String_Constant>(list_item)) {
            continue;
          }
        }
      }
      if (items_output) {
        append_string(sep);
      }
      if (items_output && sep != " ")
        append_optional_space();
      list_item->perform(this);
      items_output = true;
    }

    in_comma_array = was_comma_array;
    in_space_array = was_space_array;

    // if the list is bracketed, always include the right bracket
    if (list->is_bracketed()) {
      if (list->separator() == SASS_COMMA && list->size() == 1) {
        append_string(",");
      }
      append_string(rbracket(list));
    }
    // probably ruby sass eqivalent of element_needs_parens
    else if (output_style() == TO_SASS &&
        list->length() == 1 &&
        !list->from_selector() &&
        !Cast<List>(list->at(0)) &&
        !Cast<Selector_List>(list->at(0))
    ) {
      append_string(",");
      append_string(rbracket(list));
    }
    else if (!in_declaration && (list->separator() == SASS_HASH ||
        (list->separator() == SASS_SPACE && in_space_array) ||
        (list->separator() == SASS_COMMA && in_comma_array)
    )) {
      append_string(rbracket(list));
    }

  }

  void Inspect::operator()(Binary_Expression_Ptr expr)
  {
    expr->left()->perform(this);
    if ( in_media_block ||
         (output_style() == INSPECT) || (
          expr->op().ws_before
          && (!expr->is_interpolant())
          && (expr->is_left_interpolant() ||
              expr->is_right_interpolant())

    )) append_string(" ");
    switch (expr->optype()) {
      case Sass_OP::AND: append_string("&&"); break;
      case Sass_OP::OR:  append_string("||");  break;
      case Sass_OP::EQ:  append_string("==");  break;
      case Sass_OP::NEQ: append_string("!=");  break;
      case Sass_OP::GT:  append_string(">");   break;
      case Sass_OP::GTE: append_string(">=");  break;
      case Sass_OP::LT:  append_string("<");   break;
      case Sass_OP::LTE: append_string("<=");  break;
      case Sass_OP::ADD: append_string("+");   break;
      case Sass_OP::SUB: append_string("-");   break;
      case Sass_OP::MUL: append_string("*");   break;
      case Sass_OP::DIV: append_string("/"); break;
      case Sass_OP::MOD: append_string("%");   break;
      default: break; // shouldn't get here
    }
    if ( in_media_block ||
         (output_style() == INSPECT) || (
          expr->op().ws_after
          && (!expr->is_interpolant())
          && (expr->is_left_interpolant() ||
              expr->is_right_interpolant())
    )) append_string(" ");
    expr->right()->perform(this);
  }

  void Inspect::operator()(Unary_Expression_Ptr expr)
  {
    if (expr->optype() == Unary_Expression::PLUS)       append_string("+");
    else if (expr->optype() == Unary_Expression::SLASH) append_string("/");
    else                                                append_string("-");
    expr->operand()->perform(this);
  }

  void Inspect::operator()(Function_Call_Ptr call)
  {
    append_token(call->name(), call);
    call->arguments()->perform(this);
  }

  void Inspect::operator()(Variable_Ptr var)
  {
    append_token(var->name(), var);
  }

  void Inspect::operator()(Number_Ptr n)
  {

    // reduce units
    n->reduce();

    std::stringstream ss;
    ss.precision(opt.precision);
    ss << std::fixed << n->value();

    std::string res = ss.str();
    int s = res.length();

    // delete trailing zeros
    for(s = s - 1; s > 0; --s)
    {
        if(res[s] == '0') {
          res.erase(s, 1);
        }
        else break;
    }

    // delete trailing decimal separator
    if(res[s] == '.') res.erase(s, 1);

    // some final cosmetics
    if (res == "0.0") res = "0";
    else if (res == "") res = "0";
    else if (res == "-0") res = "0";
    else if (res == "-0.0") res = "0";
    else if (opt.output_style == COMPRESSED)
    {
      if (n->zero()) {
        // check if handling negative nr
        size_t off = res[0] == '-' ? 1 : 0;
        // remove leading zero from floating point in compressed mode
        if (res[off] == '0' && res[off+1] == '.') res.erase(off, 1);
      }
    }

    // add unit now
    res += n->unit();

    // output the final token
    append_token(res, n);
  }

  // helper function for serializing colors
  template <size_t range>
  static double cap_channel(double c) {
    if      (c > range) return range;
    else if (c < 0)     return 0;
    else                return c;
  }

  void Inspect::operator()(Color_RGBA_Ptr c)
  {
    // output the final token
    std::stringstream ss;

    // original color name
    // maybe an unknown token
    std::string name = c->disp();

    // resolved color
    std::string res_name = name;

    double r = Sass::round(cap_channel<0xff>(c->r()), opt.precision);
    double g = Sass::round(cap_channel<0xff>(c->g()), opt.precision);
    double b = Sass::round(cap_channel<0xff>(c->b()), opt.precision);
    double a = cap_channel<1>   (c->a());

    // get color from given name (if one was given at all)
    if (name != "" && name_to_color(name)) {
      Color_RGBA_Ptr_Const n = name_to_color(name);
      r = Sass::round(cap_channel<0xff>(n->r()), opt.precision);
      g = Sass::round(cap_channel<0xff>(n->g()), opt.precision);
      b = Sass::round(cap_channel<0xff>(n->b()), opt.precision);
      a = cap_channel<1>   (n->a());
    }
    // otherwise get the possible resolved color name
    else {
      double numval = r * 0x10000 + g * 0x100 + b;
      if (color_to_name(numval))
        res_name = color_to_name(numval);
    }

    std::stringstream hexlet;
    // dart sass compressed all colors in regular css always
    // ruby sass and libsass does it only when not delayed
    // since color math is going to be removed, this can go too
    bool compressed = opt.output_style == COMPRESSED;
    hexlet << '#' << std::setw(1) << std::setfill('0');
    // create a short color hexlet if there is any need for it
    if (compressed && is_color_doublet(r, g, b) && a == 1) {
      hexlet << std::hex << std::setw(1) << (static_cast<unsigned long>(r) >> 4);
      hexlet << std::hex << std::setw(1) << (static_cast<unsigned long>(g) >> 4);
      hexlet << std::hex << std::setw(1) << (static_cast<unsigned long>(b) >> 4);
    } else {
      hexlet << std::hex << std::setw(2) << static_cast<unsigned long>(r);
      hexlet << std::hex << std::setw(2) << static_cast<unsigned long>(g);
      hexlet << std::hex << std::setw(2) << static_cast<unsigned long>(b);
    }

    if (compressed && !c->is_delayed()) name = "";
    if (opt.output_style == INSPECT && a >= 1) {
      append_token(hexlet.str(), c);
      return;
    }

    // retain the originally specified color definition if unchanged
    if (name != "") {
      ss << name;
    }
    else if (a >= 1) {
      if (res_name != "") {
        if (compressed && hexlet.str().size() < res_name.size()) {
          ss << hexlet.str();
        } else {
          ss << res_name;
        }
      }
      else {
        ss << hexlet.str();
      }
    }
    else {
      ss << "rgba(";
      ss << static_cast<unsigned long>(r) << ",";
      if (!compressed) ss << " ";
      ss << static_cast<unsigned long>(g) << ",";
      if (!compressed) ss << " ";
      ss << static_cast<unsigned long>(b) << ",";
      if (!compressed) ss << " ";
      ss << a << ')';
    }

    append_token(ss.str(), c);

  }

  void Inspect::operator()(Color_HSLA_Ptr c)
  {
    Color_RGBA_Obj rgba = c->toRGBA();
    operator()(rgba);
  }

  void Inspect::operator()(Boolean_Ptr b)
  {
    // output the final token
    append_token(b->value() ? "true" : "false", b);
  }

  void Inspect::operator()(String_Schema_Ptr ss)
  {
    // Evaluation should turn these into String_Constants,
    // so this method is only for inspection purposes.
    for (size_t i = 0, L = ss->length(); i < L; ++i) {
      if ((*ss)[i]->is_interpolant()) append_string("#{");
      (*ss)[i]->perform(this);
      if ((*ss)[i]->is_interpolant()) append_string("}");
    }
  }

  void Inspect::operator()(String_Constant_Ptr s)
  {
    append_token(s->value(), s);
  }

  void Inspect::operator()(String_Quoted_Ptr s)
  {
    if (const char q = s->quote_mark()) {
      append_token(quote(s->value(), q), s);
    } else {
      append_token(s->value(), s);
    }
  }

  void Inspect::operator()(Custom_Error_Ptr e)
  {
    append_token(e->message(), e);
  }

  void Inspect::operator()(Custom_Warning_Ptr w)
  {
    append_token(w->message(), w);
  }

  void Inspect::operator()(Supports_Operator_Ptr so)
  {

    if (so->needs_parens(so->left())) append_string("(");
    so->left()->perform(this);
    if (so->needs_parens(so->left())) append_string(")");

    if (so->operand() == Supports_Operator::AND) {
      append_mandatory_space();
      append_token("and", so);
      append_mandatory_space();
    } else if (so->operand() == Supports_Operator::OR) {
      append_mandatory_space();
      append_token("or", so);
      append_mandatory_space();
    }

    if (so->needs_parens(so->right())) append_string("(");
    so->right()->perform(this);
    if (so->needs_parens(so->right())) append_string(")");
  }

  void Inspect::operator()(Supports_Negation_Ptr sn)
  {
    append_token("not", sn);
    append_mandatory_space();
    if (sn->needs_parens(sn->condition())) append_string("(");
    sn->condition()->perform(this);
    if (sn->needs_parens(sn->condition())) append_string(")");
  }

  void Inspect::operator()(Supports_Declaration_Ptr sd)
  {
    append_string("(");
    sd->feature()->perform(this);
    append_string(": ");
    sd->value()->perform(this);
    append_string(")");
  }

  void Inspect::operator()(Supports_Interpolation_Ptr sd)
  {
    sd->value()->perform(this);
  }

  void Inspect::operator()(Media_Query_Ptr mq)
  {
    size_t i = 0;
    if (mq->media_type()) {
      if      (mq->is_negated())    append_string("not ");
      else if (mq->is_restricted()) append_string("only ");
      mq->media_type()->perform(this);
    }
    else {
      (*mq)[i++]->perform(this);
    }
    for (size_t L = mq->length(); i < L; ++i) {
      append_string(" and ");
      (*mq)[i]->perform(this);
    }
  }

  void Inspect::operator()(Media_Query_Expression_Ptr mqe)
  {
    if (mqe->is_interpolated()) {
      mqe->feature()->perform(this);
    }
    else {
      append_string("(");
      mqe->feature()->perform(this);
      if (mqe->value()) {
        append_string(": "); // verified
        mqe->value()->perform(this);
      }
      append_string(")");
    }
  }

  void Inspect::operator()(At_Root_Query_Ptr ae)
  {
    if (ae->feature()) {
      append_string("(");
      ae->feature()->perform(this);
      if (ae->value()) {
        append_colon_separator();
        ae->value()->perform(this);
      }
      append_string(")");
    }
  }

  void Inspect::operator()(Function_Ptr f)
  {
    append_token("get-function", f);
    append_string("(");
    append_string(quote(f->name()));
    append_string(")");
  }

  void Inspect::operator()(Null_Ptr n)
  {
    // output the final token
    append_token("null", n);
  }

  // parameters and arguments
  void Inspect::operator()(Parameter_Ptr p)
  {
    append_token(p->name(), p);
    if (p->default_value()) {
      append_colon_separator();
      p->default_value()->perform(this);
    }
    else if (p->is_rest_parameter()) {
      append_string("...");
    }
  }

  void Inspect::operator()(Parameters_Ptr p)
  {
    append_string("(");
    if (!p->empty()) {
      (*p)[0]->perform(this);
      for (size_t i = 1, L = p->length(); i < L; ++i) {
        append_comma_separator();
        (*p)[i]->perform(this);
      }
    }
    append_string(")");
  }

  void Inspect::operator()(Argument_Ptr a)
  {
    if (!a->name().empty()) {
      append_token(a->name(), a);
      append_colon_separator();
    }
    if (!a->value()) return;
    // Special case: argument nulls can be ignored
    if (a->value()->concrete_type() == Expression::NULL_VAL) {
      return;
    }
    if (a->value()->concrete_type() == Expression::STRING) {
      String_Constant_Ptr s = Cast<String_Constant>(a->value());
      if (s) s->perform(this);
    } else {
      a->value()->perform(this);
    }
    if (a->is_rest_argument()) {
      append_string("...");
    }
  }

  void Inspect::operator()(Arguments_Ptr a)
  {
    append_string("(");
    if (!a->empty()) {
      (*a)[0]->perform(this);
      for (size_t i = 1, L = a->length(); i < L; ++i) {
        append_string(", "); // verified
        // Sass Bug? append_comma_separator();
        (*a)[i]->perform(this);
      }
    }
    append_string(")");
  }

  void Inspect::operator()(Selector_Schema_Ptr s)
  {
    s->contents()->perform(this);
  }

  void Inspect::operator()(Parent_Selector_Ptr p)
  {
    if (p->real()) append_string("&");
  }

  void Inspect::operator()(Placeholder_Selector_Ptr s)
  {
    append_token(s->name(), s);
    if (s->has_line_break()) append_optional_linefeed();
    if (s->has_line_break()) append_indentation();

  }

  void Inspect::operator()(Type_Selector_Ptr s)
  {
    append_token(s->ns_name(), s);
  }

  void Inspect::operator()(Class_Selector_Ptr s)
  {
    append_token(s->ns_name(), s);
    if (s->has_line_break()) append_optional_linefeed();
    if (s->has_line_break()) append_indentation();
  }

  void Inspect::operator()(Id_Selector_Ptr s)
  {
    append_token(s->ns_name(), s);
    if (s->has_line_break()) append_optional_linefeed();
    if (s->has_line_break()) append_indentation();
  }

  void Inspect::operator()(Attribute_Selector_Ptr s)
  {
    append_string("[");
    add_open_mapping(s);
    append_token(s->ns_name(), s);
    if (!s->matcher().empty()) {
      append_string(s->matcher());
      if (s->value() && *s->value()) {
        s->value()->perform(this);
      }
    }
    add_close_mapping(s);
    if (s->modifier() != 0) {
      append_mandatory_space();
      append_char(s->modifier());
    }
    append_string("]");
  }

  void Inspect::operator()(Pseudo_Selector_Ptr s)
  {
    append_token(s->ns_name(), s);
    if (s->expression()) {
      append_string("(");
      s->expression()->perform(this);
      append_string(")");
    }
  }

  void Inspect::operator()(Wrapped_Selector_Ptr s)
  {
    if (s->name() == " ") {
      append_string("");
    } else {
      bool was = in_wrapped;
      in_wrapped = true;
      append_token(s->name(), s);
      append_string("(");
      bool was_comma_array = in_comma_array;
      in_comma_array = false;
      s->selector()->perform(this);
      in_comma_array = was_comma_array;
      append_string(")");
      in_wrapped = was;
    }
  }

  void Inspect::operator()(Compound_Selector_Ptr s)
  {
    for (size_t i = 0, L = s->length(); i < L; ++i) {
      (*s)[i]->perform(this);
    }
    if (s->has_line_break()) {
      if (output_style() != COMPACT) {
        append_optional_linefeed();
      }
    }
  }

  void Inspect::operator()(Complex_Selector_Ptr c)
  {
    Compound_Selector_Obj      head = c->head();
    Complex_Selector_Obj            tail = c->tail();
    Complex_Selector::Combinator comb = c->combinator();

    if (comb == Complex_Selector::ANCESTOR_OF && (!head || head->empty())) {
      if (tail) tail->perform(this);
      return;
    }

    if (c->has_line_feed()) {
      if (!(c->has_parent_ref())) {
        append_optional_linefeed();
        append_indentation();
      }
    }

    if (head && head->length() != 0) head->perform(this);
    bool is_empty = !head || head->length() == 0 || head->is_empty_reference();
    bool is_tail = head && !head->is_empty_reference() && tail;
    if (output_style() == COMPRESSED && comb != Complex_Selector::ANCESTOR_OF) scheduled_space = 0;

    switch (comb) {
      case Complex_Selector::ANCESTOR_OF:
        if (is_tail) append_mandatory_space();
      break;
      case Complex_Selector::PARENT_OF:
        append_optional_space();
        append_string(">");
        append_optional_space();
      break;
      case Complex_Selector::ADJACENT_TO:
        append_optional_space();
        append_string("+");
        append_optional_space();
      break;
      case Complex_Selector::REFERENCE:
        append_mandatory_space();
        append_string("/");
        if (c->reference()) c->reference()->perform(this);
        append_string("/");
        append_mandatory_space();
      break;
      case Complex_Selector::PRECEDES:
        if (is_empty) append_optional_space();
        else append_mandatory_space();
        append_string("~");
        if (tail) append_mandatory_space();
        else append_optional_space();
      break;
      default: break;
    }
    if (tail && comb != Complex_Selector::ANCESTOR_OF) {
      if (c->has_line_break()) append_optional_linefeed();
    }
    if (tail) tail->perform(this);
    if (!tail && c->has_line_break()) {
      if (output_style() == COMPACT) {
        append_mandatory_space();
      }
    }
  }

  void Inspect::operator()(Selector_List_Ptr g)
  {

    if (g->empty()) {
      if (output_style() == TO_SASS) {
        append_token("()", g);
      }
      return;
    }


    bool was_comma_array = in_comma_array;
    // probably ruby sass eqivalent of element_needs_parens
    if (output_style() == TO_SASS && g->length() == 1 &&
      (!Cast<List>((*g)[0]) &&
       !Cast<Selector_List>((*g)[0]))) {
      append_string("(");
    }
    else if (!in_declaration && in_comma_array) {
      append_string("(");
    }

    if (in_declaration) in_comma_array = true;

    for (size_t i = 0, L = g->length(); i < L; ++i) {
      if (!in_wrapped && i == 0) append_indentation();
      if ((*g)[i] == 0) continue;
      schedule_mapping(g->at(i)->last());
      // add_open_mapping((*g)[i]->last());
      (*g)[i]->perform(this);
      // add_close_mapping((*g)[i]->last());
      if (i < L - 1) {
        scheduled_space = 0;
        append_comma_separator();
      }
    }

    in_comma_array = was_comma_array;
    // probably ruby sass eqivalent of element_needs_parens
    if (output_style() == TO_SASS && g->length() == 1 &&
      (!Cast<List>((*g)[0]) &&
       !Cast<Selector_List>((*g)[0]))) {
      append_string(",)");
    }
    else if (!in_declaration && in_comma_array) {
      append_string(")");
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.

#include <cstdlib>


namespace Sass {

  Eval::Eval(Expand& exp)
  : exp(exp),
    ctx(exp.ctx),
    traces(exp.traces),
    force(false),
    is_in_comment(false),
    is_in_selector_schema(false)
  {
    bool_true = SASS_MEMORY_NEW(Boolean, "[NA]", true);
    bool_false = SASS_MEMORY_NEW(Boolean, "[NA]", false);
  }
  Eval::~Eval() { }

  Env* Eval::environment()
  {
    return exp.environment();
  }

  const std::string Eval::cwd()
  {
    return ctx.cwd();
  }

  struct Sass_Inspect_Options& Eval::options()
  {
    return ctx.c_options;
  }

  struct Sass_Compiler* Eval::compiler()
  {
    return ctx.c_compiler;
  }

  EnvStack& Eval::env_stack()
  {
    return exp.env_stack;
  }

  Selector_List_Obj Eval::selector()
  {
    return exp.selector();
  }

  std::vector<Sass_Callee>& Eval::callee_stack()
  {
    return ctx.callee_stack;
  }


  SelectorStack& Eval::selector_stack()
  {
    return exp.selector_stack;
  }

  bool& Eval::old_at_root_without_rule()
  {
    return exp.old_at_root_without_rule;
  }


  Expression_Ptr Eval::operator()(Block_Ptr b)
  {
    Expression_Ptr val = 0;
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      val = b->at(i)->perform(this);
      if (val) return val;
    }
    return val;
  }

  Expression_Ptr Eval::operator()(Assignment_Ptr a)
  {
    Env* env = environment();
    std::string var(a->variable());
    if (a->is_global()) {
      if (a->is_default()) {
        if (env->has_global(var)) {
          Expression_Ptr e = Cast<Expression>(env->get_global(var));
          if (!e || e->concrete_type() == Expression::NULL_VAL) {
            env->set_global(var, a->value()->perform(this));
          }
        }
        else {
          env->set_global(var, a->value()->perform(this));
        }
      }
      else {
        env->set_global(var, a->value()->perform(this));
      }
    }
    else if (a->is_default()) {
      if (env->has_lexical(var)) {
        auto cur = env;
        while (cur && cur->is_lexical()) {
          if (cur->has_local(var)) {
            if (AST_Node_Obj node = cur->get_local(var)) {
              Expression_Ptr e = Cast<Expression>(node);
              if (!e || e->concrete_type() == Expression::NULL_VAL) {
                cur->set_local(var, a->value()->perform(this));
              }
            }
            else {
              throw std::runtime_error("Env not in sync");
            }
            return 0;
          }
          cur = cur->parent();
        }
        throw std::runtime_error("Env not in sync");
      }
      else if (env->has_global(var)) {
        if (AST_Node_Obj node = env->get_global(var)) {
          Expression_Ptr e = Cast<Expression>(node);
          if (!e || e->concrete_type() == Expression::NULL_VAL) {
            env->set_global(var, a->value()->perform(this));
          }
        }
      }
      else if (env->is_lexical()) {
        env->set_local(var, a->value()->perform(this));
      }
      else {
        env->set_local(var, a->value()->perform(this));
      }
    }
    else {
      env->set_lexical(var, a->value()->perform(this));
    }
    return 0;
  }

  Expression_Ptr Eval::operator()(If_Ptr i)
  {
    Expression_Obj rv;
    Env env(environment());
    env_stack().push_back(&env);
    Expression_Obj cond = i->predicate()->perform(this);
    if (!cond->is_false()) {
      rv = i->block()->perform(this);
    }
    else {
      Block_Obj alt = i->alternative();
      if (alt) rv = alt->perform(this);
    }
    env_stack().pop_back();
    return rv.detach();
  }

  // For does not create a new env scope
  // But iteration vars are reset afterwards
  Expression_Ptr Eval::operator()(For_Ptr f)
  {
    std::string variable(f->variable());
    Expression_Obj low = f->lower_bound()->perform(this);
    if (low->concrete_type() != Expression::NUMBER) {
      traces.push_back(Backtrace(low->pstate()));
      throw Exception::TypeMismatch(traces, *low, "integer");
    }
    Expression_Obj high = f->upper_bound()->perform(this);
    if (high->concrete_type() != Expression::NUMBER) {
      traces.push_back(Backtrace(high->pstate()));
      throw Exception::TypeMismatch(traces, *high, "integer");
    }
    Number_Obj sass_start = Cast<Number>(low);
    Number_Obj sass_end = Cast<Number>(high);
    // check if units are valid for sequence
    if (sass_start->unit() != sass_end->unit()) {
      std::stringstream msg; msg << "Incompatible units: '"
        << sass_end->unit() << "' and '"
        << sass_start->unit() << "'.";
      error(msg.str(), low->pstate(), traces);
    }
    double start = sass_start->value();
    double end = sass_end->value();
    // only create iterator once in this environment
    Env env(environment(), true);
    env_stack().push_back(&env);
    Block_Obj body = f->block();
    Expression_Ptr val = 0;
    if (start < end) {
      if (f->is_inclusive()) ++end;
      for (double i = start;
           i < end;
           ++i) {
        Number_Obj it = SASS_MEMORY_NEW(Number, low->pstate(), i, sass_end->unit());
        env.set_local(variable, it);
        val = body->perform(this);
        if (val) break;
      }
    } else {
      if (f->is_inclusive()) --end;
      for (double i = start;
           i > end;
           --i) {
        Number_Obj it = SASS_MEMORY_NEW(Number, low->pstate(), i, sass_end->unit());
        env.set_local(variable, it);
        val = body->perform(this);
        if (val) break;
      }
    }
    env_stack().pop_back();
    return val;
  }

  // Eval does not create a new env scope
  // But iteration vars are reset afterwards
  Expression_Ptr Eval::operator()(Each_Ptr e)
  {
    std::vector<std::string> variables(e->variables());
    Expression_Obj expr = e->list()->perform(this);
    Env env(environment(), true);
    env_stack().push_back(&env);
    List_Obj list;
    Map_Ptr map = nullptr;
    if (expr->concrete_type() == Expression::MAP) {
      map = Cast<Map>(expr);
    }
    else if (Selector_List_Ptr ls = Cast<Selector_List>(expr)) {
      Listize listize;
      Expression_Obj rv = ls->perform(&listize);
      list = Cast<List>(rv);
    }
    else if (expr->concrete_type() != Expression::LIST) {
      list = SASS_MEMORY_NEW(List, expr->pstate(), 1, SASS_COMMA);
      list->append(expr);
    }
    else {
      list = Cast<List>(expr);
    }

    Block_Obj body = e->block();
    Expression_Obj val;

    if (map) {
      for (Expression_Obj key : map->keys()) {
        Expression_Obj value = map->at(key);

        if (variables.size() == 1) {
          List_Ptr variable = SASS_MEMORY_NEW(List, map->pstate(), 2, SASS_SPACE);
          variable->append(key);
          variable->append(value);
          env.set_local(variables[0], variable);
        } else {
          env.set_local(variables[0], key);
          env.set_local(variables[1], value);
        }

        val = body->perform(this);
        if (val) break;
      }
    }
    else {
      if (list->length() == 1 && Cast<Selector_List>(list)) {
        list = Cast<List>(list);
      }
      for (size_t i = 0, L = list->length(); i < L; ++i) {
        Expression_Ptr item = list->at(i);
        // unwrap value if the expression is an argument
        if (Argument_Ptr arg = Cast<Argument>(item)) item = arg->value();
        // check if we got passed a list of args (investigate)
        if (List_Ptr scalars = Cast<List>(item)) {
          if (variables.size() == 1) {
            Expression_Ptr var = scalars;
            env.set_local(variables[0], var);
          } else {
            // XXX: this is never hit via spec tests
            for (size_t j = 0, K = variables.size(); j < K; ++j) {
              Expression_Ptr res = j >= scalars->length()
                ? SASS_MEMORY_NEW(Null, expr->pstate())
                : scalars->at(j);
              env.set_local(variables[j], res);
            }
          }
        } else {
          if (variables.size() > 0) {
            env.set_local(variables.at(0), item);
            for (size_t j = 1, K = variables.size(); j < K; ++j) {
              // XXX: this is never hit via spec tests
              Expression_Ptr res = SASS_MEMORY_NEW(Null, expr->pstate());
              env.set_local(variables[j], res);
            }
          }
        }
        val = body->perform(this);
        if (val) break;
      }
    }
    env_stack().pop_back();
    return val.detach();
  }

  Expression_Ptr Eval::operator()(While_Ptr w)
  {
    Expression_Obj pred = w->predicate();
    Block_Obj body = w->block();
    Env env(environment(), true);
    env_stack().push_back(&env);
    Expression_Obj cond = pred->perform(this);
    while (!cond->is_false()) {
      Expression_Obj val = body->perform(this);
      if (val) {
        env_stack().pop_back();
        return val.detach();
      }
      cond = pred->perform(this);
    }
    env_stack().pop_back();
    return 0;
  }

  Expression_Ptr Eval::operator()(Return_Ptr r)
  {
    return r->value()->perform(this);
  }

  Expression_Ptr Eval::operator()(Warning_Ptr w)
  {
    Sass_Output_Style outstyle = options().output_style;
    options().output_style = NESTED;
    Expression_Obj message = w->message()->perform(this);
    Env* env = environment();

    // try to use generic function
    if (env->has("@warn[f]")) {

      // add call stack entry
      callee_stack().push_back({
        "@warn",
        w->pstate().path,
        w->pstate().line + 1,
        w->pstate().column + 1,
        SASS_CALLEE_FUNCTION,
        { env }
      });

      Definition_Ptr def = Cast<Definition>((*env)["@warn[f]"]);
      // Block_Obj          body   = def->block();
      // Native_Function func   = def->native_function();
      Sass_Function_Entry c_function = def->c_function();
      Sass_Function_Fn c_func = sass_function_get_function(c_function);

      AST2C ast2c;
      union Sass_Value* c_args = sass_make_list(1, SASS_COMMA, false);
      sass_list_set_value(c_args, 0, message->perform(&ast2c));
      union Sass_Value* c_val = c_func(c_args, c_function, compiler());
      options().output_style = outstyle;
      callee_stack().pop_back();
      sass_delete_value(c_args);
      sass_delete_value(c_val);
      return 0;

    }

    std::string result(unquote(message->to_sass()));
    std::cerr << "WARNING: " << result << std::endl;
    traces.push_back(Backtrace(w->pstate()));
    std::cerr << traces_to_string(traces, "         ");
    std::cerr << std::endl;
    options().output_style = outstyle;
    traces.pop_back();
    return 0;
  }

  Expression_Ptr Eval::operator()(Error_Ptr e)
  {
    Sass_Output_Style outstyle = options().output_style;
    options().output_style = NESTED;
    Expression_Obj message = e->message()->perform(this);
    Env* env = environment();

    // try to use generic function
    if (env->has("@error[f]")) {

      // add call stack entry
      callee_stack().push_back({
        "@error",
        e->pstate().path,
        e->pstate().line + 1,
        e->pstate().column + 1,
        SASS_CALLEE_FUNCTION,
        { env }
      });

      Definition_Ptr def = Cast<Definition>((*env)["@error[f]"]);
      // Block_Obj          body   = def->block();
      // Native_Function func   = def->native_function();
      Sass_Function_Entry c_function = def->c_function();
      Sass_Function_Fn c_func = sass_function_get_function(c_function);

      AST2C ast2c;
      union Sass_Value* c_args = sass_make_list(1, SASS_COMMA, false);
      sass_list_set_value(c_args, 0, message->perform(&ast2c));
      union Sass_Value* c_val = c_func(c_args, c_function, compiler());
      options().output_style = outstyle;
      callee_stack().pop_back();
      sass_delete_value(c_args);
      sass_delete_value(c_val);
      return 0;

    }

    std::string result(unquote(message->to_sass()));
    options().output_style = outstyle;
    error(result, e->pstate(), traces);
    return 0;
  }

  Expression_Ptr Eval::operator()(Debug_Ptr d)
  {
    Sass_Output_Style outstyle = options().output_style;
    options().output_style = NESTED;
    Expression_Obj message = d->value()->perform(this);
    Env* env = environment();

    // try to use generic function
    if (env->has("@debug[f]")) {

      // add call stack entry
      callee_stack().push_back({
        "@debug",
        d->pstate().path,
        d->pstate().line + 1,
        d->pstate().column + 1,
        SASS_CALLEE_FUNCTION,
        { env }
      });

      Definition_Ptr def = Cast<Definition>((*env)["@debug[f]"]);
      // Block_Obj          body   = def->block();
      // Native_Function func   = def->native_function();
      Sass_Function_Entry c_function = def->c_function();
      Sass_Function_Fn c_func = sass_function_get_function(c_function);

      AST2C ast2c;
      union Sass_Value* c_args = sass_make_list(1, SASS_COMMA, false);
      sass_list_set_value(c_args, 0, message->perform(&ast2c));
      union Sass_Value* c_val = c_func(c_args, c_function, compiler());
      options().output_style = outstyle;
      callee_stack().pop_back();
      sass_delete_value(c_args);
      sass_delete_value(c_val);
      return 0;

    }

    std::string result(unquote(message->to_sass()));
    std::string abs_path(Sass::File::rel2abs(d->pstate().path, cwd(), cwd()));
    std::string rel_path(Sass::File::abs2rel(d->pstate().path, cwd(), cwd()));
    std::string output_path(Sass::File::path_for_console(rel_path, abs_path, d->pstate().path));
    options().output_style = outstyle;

    std::cerr << output_path << ":" << d->pstate().line+1 << " DEBUG: " << result;
    std::cerr << std::endl;
    return 0;
  }

  Expression_Ptr Eval::operator()(List_Ptr l)
  {
    // special case for unevaluated map
    if (l->separator() == SASS_HASH) {
      Map_Obj lm = SASS_MEMORY_NEW(Map,
                                l->pstate(),
                                l->length() / 2);
      for (size_t i = 0, L = l->length(); i < L; i += 2)
      {
        Expression_Obj key = (*l)[i+0]->perform(this);
        Expression_Obj val = (*l)[i+1]->perform(this);
        // make sure the color key never displays its real name
        key->is_delayed(true); // verified
        *lm << std::make_pair(key, val);
      }
      if (lm->has_duplicate_key()) {
        traces.push_back(Backtrace(l->pstate()));
        throw Exception::DuplicateKeyError(traces, *lm, *l);
      }

      lm->is_interpolant(l->is_interpolant());
      return lm->perform(this);
    }
    // check if we should expand it
    if (l->is_expanded()) return l;
    // regular case for unevaluated lists
    List_Obj ll = SASS_MEMORY_NEW(List,
                               l->pstate(),
                               l->length(),
                               l->separator(),
                               l->is_arglist(),
                               l->is_bracketed());
    for (size_t i = 0, L = l->length(); i < L; ++i) {
      ll->append((*l)[i]->perform(this));
    }
    ll->is_interpolant(l->is_interpolant());
    ll->from_selector(l->from_selector());
    ll->is_expanded(true);
    return ll.detach();
  }

  Expression_Ptr Eval::operator()(Map_Ptr m)
  {
    if (m->is_expanded()) return m;

    // make sure we're not starting with duplicate keys.
    // the duplicate key state will have been set in the parser phase.
    if (m->has_duplicate_key()) {
      traces.push_back(Backtrace(m->pstate()));
      throw Exception::DuplicateKeyError(traces, *m, *m);
    }

    Map_Obj mm = SASS_MEMORY_NEW(Map,
                                m->pstate(),
                                m->length());
    for (auto key : m->keys()) {
      Expression_Ptr ex_key = key->perform(this);
      Expression_Ptr ex_val = m->at(key);
      if (ex_val == NULL) continue;
      ex_val = ex_val->perform(this);
      *mm << std::make_pair(ex_key, ex_val);
    }

    // check the evaluated keys aren't duplicates.
    if (mm->has_duplicate_key()) {
      traces.push_back(Backtrace(m->pstate()));
      throw Exception::DuplicateKeyError(traces, *mm, *m);
    }

    mm->is_expanded(true);
    return mm.detach();
  }

  Expression_Ptr Eval::operator()(Binary_Expression_Ptr b_in)
  {

    Expression_Obj lhs = b_in->left();
    Expression_Obj rhs = b_in->right();
    enum Sass_OP op_type = b_in->optype();

    if (op_type == Sass_OP::AND) {
      // LOCAL_FLAG(force, true);
      lhs = lhs->perform(this);
      if (!*lhs) return lhs.detach();
      return rhs->perform(this);
    }
    else if (op_type == Sass_OP::OR) {
      // LOCAL_FLAG(force, true);
      lhs = lhs->perform(this);
      if (*lhs) return lhs.detach();
      return rhs->perform(this);
    }

    // Evaluate variables as early o
    while (Variable_Ptr l_v = Cast<Variable>(lhs)) {
      lhs = operator()(l_v);
    }
    while (Variable_Ptr r_v = Cast<Variable>(rhs)) {
      rhs = operator()(r_v);
    }

    Binary_Expression_Obj b = b_in;

    // Evaluate sub-expressions early on
    while (Binary_Expression_Ptr l_b = Cast<Binary_Expression>(lhs)) {
      if (!force && l_b->is_delayed()) break;
      lhs = operator()(l_b);
    }
    while (Binary_Expression_Ptr r_b = Cast<Binary_Expression>(rhs)) {
      if (!force && r_b->is_delayed()) break;
      rhs = operator()(r_b);
    }

    // don't eval delayed expressions (the '/' when used as a separator)
    if (!force && op_type == Sass_OP::DIV && b->is_delayed()) {
      b->right(b->right()->perform(this));
      b->left(b->left()->perform(this));
      return b.detach();
    }

    // specific types we know are final
    // handle them early to avoid overhead
    if (Number_Ptr l_n = Cast<Number>(lhs)) {
      // lhs is number and rhs is number
      if (Number_Ptr r_n = Cast<Number>(rhs)) {
        try {
          switch (op_type) {
            case Sass_OP::EQ: return *l_n == *r_n ? bool_true : bool_false;
            case Sass_OP::NEQ: return *l_n == *r_n ? bool_false : bool_true;
            case Sass_OP::LT: return *l_n < *r_n ? bool_true : bool_false;
            case Sass_OP::GTE: return *l_n < *r_n ? bool_false : bool_true;
            case Sass_OP::LTE: return *l_n < *r_n || *l_n == *r_n ? bool_true : bool_false;
            case Sass_OP::GT: return *l_n < *r_n || *l_n == *r_n ? bool_false : bool_true;
            case Sass_OP::ADD: case Sass_OP::SUB: case Sass_OP::MUL: case Sass_OP::DIV: case Sass_OP::MOD:
              return Operators::op_numbers(op_type, *l_n, *r_n, options(), b_in->pstate());
            default: break;
          }
        }
        catch (Exception::OperationError& err)
        {
          traces.push_back(Backtrace(b_in->pstate()));
          throw Exception::SassValueError(traces, b_in->pstate(), err);
        }
      }
      // lhs is number and rhs is color
      // Todo: allow to work with HSLA colors
      else if (Color_Ptr r_col = Cast<Color>(rhs)) {
        Color_RGBA_Obj r_c = r_col->toRGBA();
        try {
          switch (op_type) {
            case Sass_OP::EQ: return *l_n == *r_c ? bool_true : bool_false;
            case Sass_OP::NEQ: return *l_n == *r_c ? bool_false : bool_true;
            case Sass_OP::ADD: case Sass_OP::SUB: case Sass_OP::MUL: case Sass_OP::DIV: case Sass_OP::MOD:
              return Operators::op_number_color(op_type, *l_n, *r_c, options(), b_in->pstate());
            default: break;
          }
        }
        catch (Exception::OperationError& err)
        {
          traces.push_back(Backtrace(b_in->pstate()));
          throw Exception::SassValueError(traces, b_in->pstate(), err);
        }
      }
    }
    else if (Color_Ptr l_col = Cast<Color>(lhs)) {
      Color_RGBA_Obj l_c = l_col->toRGBA();
      // lhs is color and rhs is color
      if (Color_Ptr r_col = Cast<Color>(rhs)) {
        Color_RGBA_Obj r_c = r_col->toRGBA();
        try {
          switch (op_type) {
            case Sass_OP::EQ: return *l_c == *r_c ? bool_true : bool_false;
            case Sass_OP::NEQ: return *l_c == *r_c ? bool_false : bool_true;
            case Sass_OP::LT: return *l_c < *r_c ? bool_true : bool_false;
            case Sass_OP::GTE: return *l_c < *r_c ? bool_false : bool_true;
            case Sass_OP::LTE: return *l_c < *r_c || *l_c == *r_c ? bool_true : bool_false;
            case Sass_OP::GT: return *l_c < *r_c || *l_c == *r_c ? bool_false : bool_true;
            case Sass_OP::ADD: case Sass_OP::SUB: case Sass_OP::MUL: case Sass_OP::DIV: case Sass_OP::MOD:
              return Operators::op_colors(op_type, *l_c, *r_c, options(), b_in->pstate());
            default: break;
          }
        }
        catch (Exception::OperationError& err)
        {
          traces.push_back(Backtrace(b_in->pstate()));
          throw Exception::SassValueError(traces, b_in->pstate(), err);
        }
      }
      // lhs is color and rhs is number
      else if (Number_Ptr r_n = Cast<Number>(rhs)) {
        try {
          switch (op_type) {
            case Sass_OP::EQ: return *l_c == *r_n ? bool_true : bool_false;
            case Sass_OP::NEQ: return *l_c == *r_n ? bool_false : bool_true;
            case Sass_OP::ADD: case Sass_OP::SUB: case Sass_OP::MUL: case Sass_OP::DIV: case Sass_OP::MOD:
              return Operators::op_color_number(op_type, *l_c, *r_n, options(), b_in->pstate());
            default: break;
          }
        }
        catch (Exception::OperationError& err)
        {
          traces.push_back(Backtrace(b_in->pstate()));
          throw Exception::SassValueError(traces, b_in->pstate(), err);
        }
      }
    }

    String_Schema_Obj ret_schema;

    // only the last item will be used to eval the binary expression
    if (String_Schema_Ptr s_l = Cast<String_Schema>(b->left())) {
      if (!s_l->has_interpolant() && (!s_l->is_right_interpolant())) {
        ret_schema = SASS_MEMORY_NEW(String_Schema, b->pstate());
        Binary_Expression_Obj bin_ex = SASS_MEMORY_NEW(Binary_Expression, b->pstate(),
                                                    b->op(), s_l->last(), b->right());
        bin_ex->is_delayed(b->left()->is_delayed() || b->right()->is_delayed()); // unverified
        for (size_t i = 0; i < s_l->length() - 1; ++i) {
          ret_schema->append(Cast<PreValue>(s_l->at(i)->perform(this)));
        }
        ret_schema->append(Cast<PreValue>(bin_ex->perform(this)));
        return ret_schema->perform(this);
      }
    }
    if (String_Schema_Ptr s_r = Cast<String_Schema>(b->right())) {

      if (!s_r->has_interpolant() && (!s_r->is_left_interpolant() || op_type == Sass_OP::DIV)) {
        ret_schema = SASS_MEMORY_NEW(String_Schema, b->pstate());
        Binary_Expression_Obj bin_ex = SASS_MEMORY_NEW(Binary_Expression, b->pstate(),
                                                    b->op(), b->left(), s_r->first());
        bin_ex->is_delayed(b->left()->is_delayed() || b->right()->is_delayed()); // verified
        ret_schema->append(Cast<PreValue>(bin_ex->perform(this)));
        for (size_t i = 1; i < s_r->length(); ++i) {
          ret_schema->append(Cast<PreValue>(s_r->at(i)->perform(this)));
        }
        return ret_schema->perform(this);
      }
    }

    // fully evaluate their values
    if (op_type == Sass_OP::EQ ||
        op_type == Sass_OP::NEQ ||
        op_type == Sass_OP::GT ||
        op_type == Sass_OP::GTE ||
        op_type == Sass_OP::LT ||
        op_type == Sass_OP::LTE)
    {
      LOCAL_FLAG(force, true);
      lhs->is_expanded(false);
      lhs->set_delayed(false);
      lhs = lhs->perform(this);
      rhs->is_expanded(false);
      rhs->set_delayed(false);
      rhs = rhs->perform(this);
    }
    else {
      lhs = lhs->perform(this);
    }

    // not a logical connective, so go ahead and eval the rhs
    rhs = rhs->perform(this);
    AST_Node_Obj lu = lhs;
    AST_Node_Obj ru = rhs;

    Expression::Type l_type;
    Expression::Type r_type;

    // Is one of the operands an interpolant?
    String_Schema_Obj s1 = Cast<String_Schema>(b->left());
    String_Schema_Obj s2 = Cast<String_Schema>(b->right());
    Binary_Expression_Obj b1 = Cast<Binary_Expression>(b->left());
    Binary_Expression_Obj b2 = Cast<Binary_Expression>(b->right());

    bool schema_op = false;

    bool force_delay = (s2 && s2->is_left_interpolant()) ||
                       (s1 && s1->is_right_interpolant()) ||
                       (b1 && b1->is_right_interpolant()) ||
                       (b2 && b2->is_left_interpolant());

    if ((s1 && s1->has_interpolants()) || (s2 && s2->has_interpolants()) || force_delay)
    {
      if (op_type == Sass_OP::DIV || op_type == Sass_OP::MUL || op_type == Sass_OP::MOD || op_type == Sass_OP::ADD || op_type == Sass_OP::SUB ||
          op_type == Sass_OP::EQ) {
        // If possible upgrade LHS to a number (for number to string compare)
        if (String_Constant_Ptr str = Cast<String_Constant>(lhs)) {
          std::string value(str->value());
          const char* start = value.c_str();
          if (Prelexer::sequence < Prelexer::dimension, Prelexer::end_of_file >(start) != 0) {
            lhs = Parser::lexed_dimension(b->pstate(), str->value());
          }
        }
        // If possible upgrade RHS to a number (for string to number compare)
        if (String_Constant_Ptr str = Cast<String_Constant>(rhs)) {
          std::string value(str->value());
          const char* start = value.c_str();
          if (Prelexer::sequence < Prelexer::dimension, Prelexer::number >(start) != 0) {
            rhs = Parser::lexed_dimension(b->pstate(), str->value());
          }
        }
      }

      To_Value to_value(ctx);
      Value_Obj v_l = Cast<Value>(lhs->perform(&to_value));
      Value_Obj v_r = Cast<Value>(rhs->perform(&to_value));

      if (force_delay) {
        std::string str("");
        str += v_l->to_string(options());
        if (b->op().ws_before) str += " ";
        str += b->separator();
        if (b->op().ws_after) str += " ";
        str += v_r->to_string(options());
        String_Constant_Ptr val = SASS_MEMORY_NEW(String_Constant, b->pstate(), str);
        val->is_interpolant(b->left()->has_interpolant());
        return val;
      }
    }

    // see if it's a relational expression
    try {
      switch(op_type) {
        case Sass_OP::EQ:  return SASS_MEMORY_NEW(Boolean, b->pstate(), Operators::eq(lhs, rhs));
        case Sass_OP::NEQ: return SASS_MEMORY_NEW(Boolean, b->pstate(), Operators::neq(lhs, rhs));
        case Sass_OP::GT:  return SASS_MEMORY_NEW(Boolean, b->pstate(), Operators::gt(lhs, rhs));
        case Sass_OP::GTE: return SASS_MEMORY_NEW(Boolean, b->pstate(), Operators::gte(lhs, rhs));
        case Sass_OP::LT:  return SASS_MEMORY_NEW(Boolean, b->pstate(), Operators::lt(lhs, rhs));
        case Sass_OP::LTE: return SASS_MEMORY_NEW(Boolean, b->pstate(), Operators::lte(lhs, rhs));
        default: break;
      }
    }
    catch (Exception::OperationError& err)
    {
      traces.push_back(Backtrace(b->pstate()));
      throw Exception::SassValueError(traces, b->pstate(), err);
    }

    l_type = lhs->concrete_type();
    r_type = rhs->concrete_type();

    // ToDo: throw error in op functions
    // ToDo: then catch and re-throw them
    Expression_Obj rv;
    try {
      ParserState pstate(b->pstate());
      if (l_type == Expression::NUMBER && r_type == Expression::NUMBER) {
        Number_Ptr l_n = Cast<Number>(lhs);
        Number_Ptr r_n = Cast<Number>(rhs);
        l_n->reduce(); r_n->reduce();
        rv = Operators::op_numbers(op_type, *l_n, *r_n, options(), pstate);
      }
      else if (l_type == Expression::NUMBER && r_type == Expression::COLOR) {
        Number_Ptr l_n = Cast<Number>(lhs);
        Color_RGBA_Obj r_c = Cast<Color>(rhs)->toRGBA();
        rv = Operators::op_number_color(op_type, *l_n, *r_c, options(), pstate);
      }
      else if (l_type == Expression::COLOR && r_type == Expression::NUMBER) {
        Color_RGBA_Obj l_c = Cast<Color>(lhs)->toRGBA();
        Number_Ptr r_n = Cast<Number>(rhs);
        rv = Operators::op_color_number(op_type, *l_c, *r_n, options(), pstate);
      }
      else if (l_type == Expression::COLOR && r_type == Expression::COLOR) {
        Color_RGBA_Obj l_c = Cast<Color>(lhs)->toRGBA();
        Color_RGBA_Obj r_c = Cast<Color>(rhs)->toRGBA();
        rv = Operators::op_colors(op_type, *l_c, *r_c, options(), pstate);
      }
      else {
        To_Value to_value(ctx);
        // this will leak if perform does not return a value!
        Value_Obj v_l = Cast<Value>(lhs->perform(&to_value));
        Value_Obj v_r = Cast<Value>(rhs->perform(&to_value));
        bool interpolant = b->is_right_interpolant() ||
                           b->is_left_interpolant() ||
                           b->is_interpolant();
        if (op_type == Sass_OP::SUB) interpolant = false;
        // if (op_type == Sass_OP::DIV) interpolant = true;
        // check for type violations
        if (l_type == Expression::MAP || l_type == Expression::FUNCTION_VAL) {
          traces.push_back(Backtrace(v_l->pstate()));
          throw Exception::InvalidValue(traces, *v_l);
        }
        if (r_type == Expression::MAP || l_type == Expression::FUNCTION_VAL) {
          traces.push_back(Backtrace(v_r->pstate()));
          throw Exception::InvalidValue(traces, *v_r);
        }
        Value_Ptr ex = Operators::op_strings(b->op(), *v_l, *v_r, options(), pstate, !interpolant); // pass true to compress
        if (String_Constant_Ptr str = Cast<String_Constant>(ex))
        {
          if (str->concrete_type() == Expression::STRING)
          {
            String_Constant_Ptr lstr = Cast<String_Constant>(lhs);
            String_Constant_Ptr rstr = Cast<String_Constant>(rhs);
            if (op_type != Sass_OP::SUB) {
              if (String_Constant_Ptr org = lstr ? lstr : rstr)
              { str->quote_mark(org->quote_mark()); }
            }
          }
        }
        ex->is_interpolant(b->is_interpolant());
        rv = ex;
      }
    }
    catch (Exception::OperationError& err)
    {
      traces.push_back(Backtrace(b->pstate()));
      // throw Exception::Base(b->pstate(), err.what());
      throw Exception::SassValueError(traces, b->pstate(), err);
    }

    if (rv) {
      if (schema_op) {
        // XXX: this is never hit via spec tests
        (*s2)[0] = rv;
        rv = s2->perform(this);
      }
    }

    return rv.detach();

  }

  Expression_Ptr Eval::operator()(Unary_Expression_Ptr u)
  {
    Expression_Obj operand = u->operand()->perform(this);
    if (u->optype() == Unary_Expression::NOT) {
      Boolean_Ptr result = SASS_MEMORY_NEW(Boolean, u->pstate(), (bool)*operand);
      result->value(!result->value());
      return result;
    }
    else if (Number_Obj nr = Cast<Number>(operand)) {
      // negate value for minus unary expression
      if (u->optype() == Unary_Expression::MINUS) {
        Number_Obj cpy = SASS_MEMORY_COPY(nr);
        cpy->value( - cpy->value() ); // negate value
        return cpy.detach(); // return the copy
      }
      else if (u->optype() == Unary_Expression::SLASH) {
        std::string str = '/' + nr->to_string(options());
        return SASS_MEMORY_NEW(String_Constant, u->pstate(), str);
      }
      // nothing for positive
      return nr.detach();
    }
    else {
      // Special cases: +/- variables which evaluate to null ouput just +/-,
      // but +/- null itself outputs the string
      if (operand->concrete_type() == Expression::NULL_VAL && Cast<Variable>(u->operand())) {
        u->operand(SASS_MEMORY_NEW(String_Quoted, u->pstate(), ""));
      }
      // Never apply unary opertions on colors @see #2140
      else if (Color_Ptr color = Cast<Color>(operand)) {
        // Use the color name if this was eval with one
        if (color->disp().length() > 0) {
          operand = SASS_MEMORY_NEW(String_Constant, operand->pstate(), color->disp());
          u->operand(operand);
        }
      }
      else {
        u->operand(operand);
      }

      return SASS_MEMORY_NEW(String_Quoted,
                             u->pstate(),
                             u->inspect());
    }
    // unreachable
    return u;
  }

  Expression_Ptr Eval::operator()(Function_Call_Ptr c)
  {
    if (traces.size() > Constants::MaxCallStack) {
        // XXX: this is never hit via spec tests
        std::ostringstream stm;
        stm << "Stack depth exceeded max of " << Constants::MaxCallStack;
        error(stm.str(), c->pstate(), traces);
    }

    if (Cast<String_Schema>(c->sname())) {
      Expression_Obj evaluated_name = c->sname()->perform(this);
      Expression_Obj evaluated_args = c->arguments()->perform(this);
      std::string str(evaluated_name->to_string());
      str += evaluated_args->to_string();
      return SASS_MEMORY_NEW(String_Constant, c->pstate(), str);
    }

    std::string name(Util::normalize_underscores(c->name()));
    std::string full_name(name + "[f]");

    // we make a clone here, need to implement that further
    Arguments_Obj args = c->arguments();

    Env* env = environment();
    if (!env->has(full_name) || (!c->via_call() && Prelexer::re_special_fun(name.c_str()))) {
      if (!env->has("*[f]")) {
        for (Argument_Obj arg : args->elements()) {
          if (List_Obj ls = Cast<List>(arg->value())) {
            if (ls->size() == 0) error("() isn't a valid CSS value.", c->pstate(), traces);
          }
        }
        args = Cast<Arguments>(args->perform(this));
        Function_Call_Obj lit = SASS_MEMORY_NEW(Function_Call,
                                             c->pstate(),
                                             c->name(),
                                             args);
        if (args->has_named_arguments()) {
          error("Function " + c->name() + " doesn't support keyword arguments", c->pstate(), traces);
        }
        String_Quoted_Ptr str = SASS_MEMORY_NEW(String_Quoted,
                                             c->pstate(),
                                             lit->to_string(options()));
        str->is_interpolant(c->is_interpolant());
        return str;
      } else {
        // call generic function
        full_name = "*[f]";
      }
    }

    // further delay for calls
    if (full_name != "call[f]") {
      args->set_delayed(false); // verified
    }
    if (full_name != "if[f]") {
      args = Cast<Arguments>(args->perform(this));
    }
    Definition_Ptr def = Cast<Definition>((*env)[full_name]);

    if (c->func()) def = c->func()->definition();

    if (def->is_overload_stub()) {
      std::stringstream ss;
      size_t L = args->length();
      // account for rest arguments
      if (args->has_rest_argument() && args->length() > 0) {
        // get the rest arguments list
        List_Ptr rest = Cast<List>(args->last()->value());
        // arguments before rest argument plus rest
        if (rest) L += rest->length() - 1;
      }
      ss << full_name << L;
      full_name = ss.str();
      std::string resolved_name(full_name);
      if (!env->has(resolved_name)) error("overloaded function `" + std::string(c->name()) + "` given wrong number of arguments", c->pstate(), traces);
      def = Cast<Definition>((*env)[resolved_name]);
    }

    Expression_Obj     result = c;
    Block_Obj          body   = def->block();
    Native_Function func   = def->native_function();
    Sass_Function_Entry c_function = def->c_function();

    if (c->is_css()) return result.detach();

    Parameters_Obj params = def->parameters();
    Env fn_env(def->environment());
    env_stack().push_back(&fn_env);

    if (func || body) {
      bind(std::string("Function"), c->name(), params, args, &fn_env, this, traces);
      std::string msg(", in function `" + c->name() + "`");
      traces.push_back(Backtrace(c->pstate(), msg));
      callee_stack().push_back({
        c->name().c_str(),
        c->pstate().path,
        c->pstate().line + 1,
        c->pstate().column + 1,
        SASS_CALLEE_FUNCTION,
        { env }
      });

      // eval the body if user-defined or special, invoke underlying CPP function if native
      if (body /* && !Prelexer::re_special_fun(name.c_str()) */) {
        result = body->perform(this);
      }
      else if (func) {
        result = func(fn_env, *env, ctx, def->signature(), c->pstate(), traces, exp.selector_stack);
      }
      if (!result) {
        error(std::string("Function ") + c->name() + " finished without @return", c->pstate(), traces);
      }
      callee_stack().pop_back();
      traces.pop_back();
    }

    // else if it's a user-defined c function
    // convert call into C-API compatible form
    else if (c_function) {
      Sass_Function_Fn c_func = sass_function_get_function(c_function);
      if (full_name == "*[f]") {
        String_Quoted_Obj str = SASS_MEMORY_NEW(String_Quoted, c->pstate(), c->name());
        Arguments_Obj new_args = SASS_MEMORY_NEW(Arguments, c->pstate());
        new_args->append(SASS_MEMORY_NEW(Argument, c->pstate(), str));
        new_args->concat(args);
        args = new_args;
      }

      // populates env with default values for params
      std::string ff(c->name());
      bind(std::string("Function"), c->name(), params, args, &fn_env, this, traces);
      std::string msg(", in function `" + c->name() + "`");
      traces.push_back(Backtrace(c->pstate(), msg));
      callee_stack().push_back({
        c->name().c_str(),
        c->pstate().path,
        c->pstate().line + 1,
        c->pstate().column + 1,
        SASS_CALLEE_C_FUNCTION,
        { env }
      });

      AST2C ast2c;
      union Sass_Value* c_args = sass_make_list(params->length(), SASS_COMMA, false);
      for(size_t i = 0; i < params->length(); i++) {
        Parameter_Obj param = params->at(i);
        std::string key = param->name();
        AST_Node_Obj node = fn_env.get_local(key);
        Expression_Obj arg = Cast<Expression>(node);
        sass_list_set_value(c_args, i, arg->perform(&ast2c));
      }
      union Sass_Value* c_val = c_func(c_args, c_function, compiler());
      if (sass_value_get_tag(c_val) == SASS_ERROR) {
        std::string message("error in C function " + c->name() + ": " + sass_error_get_message(c_val));
        sass_delete_value(c_val);
        sass_delete_value(c_args);
        error(message, c->pstate(), traces);
      } else if (sass_value_get_tag(c_val) == SASS_WARNING) {
        std::string message("warning in C function " + c->name() + ": " + sass_warning_get_message(c_val));
        sass_delete_value(c_val);
        sass_delete_value(c_args);
        error(message, c->pstate(), traces);
      }
      result = c2ast(c_val, traces, c->pstate());

      callee_stack().pop_back();
      traces.pop_back();
      sass_delete_value(c_args);
      if (c_val != c_args)
        sass_delete_value(c_val);
    }

    // link back to function definition
    // only do this for custom functions
    if (result->pstate().file == std::string::npos)
      result->pstate(c->pstate());

    result = result->perform(this);
    result->is_interpolant(c->is_interpolant());
    env_stack().pop_back();
    return result.detach();
  }

  Expression_Ptr Eval::operator()(Variable_Ptr v)
  {
    Expression_Obj value;
    Env* env = environment();
    const std::string& name(v->name());
    EnvResult rv(env->find(name));
    if (rv.found) value = static_cast<Expression*>(rv.it->second.ptr());
    else error("Undefined variable: \"" + v->name() + "\".", v->pstate(), traces);
    if (Argument_Ptr arg = Cast<Argument>(value)) value = arg->value();
    if (Number_Ptr nr = Cast<Number>(value)) nr->zero(true); // force flag
    value->is_interpolant(v->is_interpolant());
    if (force) value->is_expanded(false);
    value->set_delayed(false); // verified
    value = value->perform(this);
    if(!force) rv.it->second = value;
    return value.detach();
  }

  Expression_Ptr Eval::operator()(Color_RGBA_Ptr c)
  {
    return c;
  }

  Expression_Ptr Eval::operator()(Color_HSLA_Ptr c)
  {
    return c;
  }

  Expression_Ptr Eval::operator()(Number_Ptr n)
  {
    return n;
  }

  Expression_Ptr Eval::operator()(Boolean_Ptr b)
  {
    return b;
  }

  void Eval::interpolation(Context& ctx, std::string& res, Expression_Obj ex, bool into_quotes, bool was_itpl) {

    bool needs_closing_brace = false;

    if (Arguments_Ptr args = Cast<Arguments>(ex)) {
      List_Ptr ll = SASS_MEMORY_NEW(List, args->pstate(), 0, SASS_COMMA);
      for(auto arg : args->elements()) {
        ll->append(arg->value());
      }
      ll->is_interpolant(args->is_interpolant());
      needs_closing_brace = true;
      res += "(";
      ex = ll;
    }
    if (Number_Ptr nr = Cast<Number>(ex)) {
      Number reduced(nr);
      reduced.reduce();
      if (!reduced.is_valid_css_unit()) {
        traces.push_back(Backtrace(nr->pstate()));
        throw Exception::InvalidValue(traces, *nr);
      }
    }
    if (Argument_Ptr arg = Cast<Argument>(ex)) {
      ex = arg->value();
    }
    if (String_Quoted_Ptr sq = Cast<String_Quoted>(ex)) {
      if (was_itpl) {
        bool was_interpolant = ex->is_interpolant();
        ex = SASS_MEMORY_NEW(String_Constant, sq->pstate(), sq->value());
        ex->is_interpolant(was_interpolant);
      }
    }

    if (Cast<Null>(ex)) { return; }

    // parent selector needs another go
    if (Cast<Parent_Selector>(ex)) {
      // XXX: this is never hit via spec tests
      ex = ex->perform(this);
    }
    // parent selector needs another go
    if (Cast<Parent_Reference>(ex)) {
      // XXX: this is never hit via spec tests
      ex = ex->perform(this);
    }

    if (List_Ptr l = Cast<List>(ex)) {
      List_Obj ll = SASS_MEMORY_NEW(List, l->pstate(), 0, l->separator());
      // this fixes an issue with bourbon sample, not really sure why
      // if (l->size() && Cast<Null>((*l)[0])) { res += ""; }
      for(Expression_Obj item : *l) {
        item->is_interpolant(l->is_interpolant());
        std::string rl(""); interpolation(ctx, rl, item, into_quotes, l->is_interpolant());
        bool is_null = Cast<Null>(item) != 0; // rl != ""
        if (!is_null) ll->append(SASS_MEMORY_NEW(String_Quoted, item->pstate(), rl));
      }
      // Check indicates that we probably should not get a list
      // here. Normally single list items are already unwrapped.
      if (l->size() > 1) {
        // string_to_output would fail "#{'_\a' '_\a'}";
        std::string str(ll->to_string(options()));
        str = read_hex_escapes(str); // read escapes
        newline_to_space(str); // replace directly
        res += str; // append to result string
      } else {
        res += (ll->to_string(options()));
      }
      ll->is_interpolant(l->is_interpolant());
    }

    // Value
    // Function_Call
    // Selector_List
    // String_Quoted
    // String_Constant
    // Parent_Selector
    // Binary_Expression
    else {
      // ex = ex->perform(this);
      if (into_quotes && ex->is_interpolant()) {
        res += evacuate_escapes(ex ? ex->to_string(options()) : "");
      } else {
        std::string str(ex ? ex->to_string(options()) : "");
        if (into_quotes) str = read_hex_escapes(str);
        res += str; // append to result string
      }
    }

    if (needs_closing_brace) res += ")";

  }

  Expression_Ptr Eval::operator()(String_Schema_Ptr s)
  {
    size_t L = s->length();
    bool into_quotes = false;
    if (L > 1) {
      if (!Cast<String_Quoted>((*s)[0]) && !Cast<String_Quoted>((*s)[L - 1])) {
      if (String_Constant_Ptr l = Cast<String_Constant>((*s)[0])) {
        if (String_Constant_Ptr r = Cast<String_Constant>((*s)[L - 1])) {
          if (r->value().size() > 0) {
            if (l->value()[0] == '"' && r->value()[r->value().size() - 1] == '"') into_quotes = true;
            if (l->value()[0] == '\'' && r->value()[r->value().size() - 1] == '\'') into_quotes = true;
          }
        }
      }
      }
    }
    bool was_quoted = false;
    bool was_interpolant = false;
    std::string res("");
    for (size_t i = 0; i < L; ++i) {
      bool is_quoted = Cast<String_Quoted>((*s)[i]) != NULL;
      if (was_quoted && !(*s)[i]->is_interpolant() && !was_interpolant) { res += " "; }
      else if (i > 0 && is_quoted && !(*s)[i]->is_interpolant() && !was_interpolant) { res += " "; }
      Expression_Obj ex = (*s)[i]->perform(this);
      interpolation(ctx, res, ex, into_quotes, ex->is_interpolant());
      was_quoted = Cast<String_Quoted>((*s)[i]) != NULL;
      was_interpolant = (*s)[i]->is_interpolant();

    }
    if (!s->is_interpolant()) {
      if (s->length() > 1 && res == "") return SASS_MEMORY_NEW(Null, s->pstate());
      return SASS_MEMORY_NEW(String_Constant, s->pstate(), res, s->css());
    }
    // string schema seems to have a special unquoting behavior (also handles "nested" quotes)
    String_Quoted_Obj str = SASS_MEMORY_NEW(String_Quoted, s->pstate(), res, 0, false, false, false, s->css());
    // if (s->is_interpolant()) str->quote_mark(0);
    // String_Constant_Ptr str = SASS_MEMORY_NEW(String_Constant, s->pstate(), res);
    if (str->quote_mark()) str->quote_mark('*');
    else if (!is_in_comment) str->value(string_to_output(str->value()));
    str->is_interpolant(s->is_interpolant());
    return str.detach();
  }


  Expression_Ptr Eval::operator()(String_Constant_Ptr s)
  {
    return s;
  }

  Expression_Ptr Eval::operator()(String_Quoted_Ptr s)
  {
    String_Quoted_Ptr str = SASS_MEMORY_NEW(String_Quoted, s->pstate(), "");
    str->value(s->value());
    str->quote_mark(s->quote_mark());
    str->is_interpolant(s->is_interpolant());
    return str;
  }

  Expression_Ptr Eval::operator()(Supports_Operator_Ptr c)
  {
    Expression_Ptr left = c->left()->perform(this);
    Expression_Ptr right = c->right()->perform(this);
    Supports_Operator_Ptr cc = SASS_MEMORY_NEW(Supports_Operator,
                                 c->pstate(),
                                 Cast<Supports_Condition>(left),
                                 Cast<Supports_Condition>(right),
                                 c->operand());
    return cc;
  }

  Expression_Ptr Eval::operator()(Supports_Negation_Ptr c)
  {
    Expression_Ptr condition = c->condition()->perform(this);
    Supports_Negation_Ptr cc = SASS_MEMORY_NEW(Supports_Negation,
                                 c->pstate(),
                                 Cast<Supports_Condition>(condition));
    return cc;
  }

  Expression_Ptr Eval::operator()(Supports_Declaration_Ptr c)
  {
    Expression_Ptr feature = c->feature()->perform(this);
    Expression_Ptr value = c->value()->perform(this);
    Supports_Declaration_Ptr cc = SASS_MEMORY_NEW(Supports_Declaration,
                              c->pstate(),
                              feature,
                              value);
    return cc;
  }

  Expression_Ptr Eval::operator()(Supports_Interpolation_Ptr c)
  {
    Expression_Ptr value = c->value()->perform(this);
    Supports_Interpolation_Ptr cc = SASS_MEMORY_NEW(Supports_Interpolation,
                            c->pstate(),
                            value);
    return cc;
  }

  Expression_Ptr Eval::operator()(At_Root_Query_Ptr e)
  {
    Expression_Obj feature = e->feature();
    feature = (feature ? feature->perform(this) : 0);
    Expression_Obj value = e->value();
    value = (value ? value->perform(this) : 0);
    Expression_Ptr ee = SASS_MEMORY_NEW(At_Root_Query,
                                     e->pstate(),
                                     Cast<String>(feature),
                                     value);
    return ee;
  }

  Media_Query_Ptr Eval::operator()(Media_Query_Ptr q)
  {
    String_Obj t = q->media_type();
    t = static_cast<String_Ptr>(t.isNull() ? 0 : t->perform(this));
    Media_Query_Obj qq = SASS_MEMORY_NEW(Media_Query,
                                      q->pstate(),
                                      t,
                                      q->length(),
                                      q->is_negated(),
                                      q->is_restricted());
    for (size_t i = 0, L = q->length(); i < L; ++i) {
      qq->append(static_cast<Media_Query_Expression_Ptr>((*q)[i]->perform(this)));
    }
    return qq.detach();
  }

  Expression_Ptr Eval::operator()(Media_Query_Expression_Ptr e)
  {
    Expression_Obj feature = e->feature();
    feature = (feature ? feature->perform(this) : 0);
    if (feature && Cast<String_Quoted>(feature)) {
      feature = SASS_MEMORY_NEW(String_Quoted,
                                  feature->pstate(),
                                  Cast<String_Quoted>(feature)->value());
    }
    Expression_Obj value = e->value();
    value = (value ? value->perform(this) : 0);
    if (value && Cast<String_Quoted>(value)) {
      // XXX: this is never hit via spec tests
      value = SASS_MEMORY_NEW(String_Quoted,
                                value->pstate(),
                                Cast<String_Quoted>(value)->value());
    }
    return SASS_MEMORY_NEW(Media_Query_Expression,
                           e->pstate(),
                           feature,
                           value,
                           e->is_interpolated());
  }

  Expression_Ptr Eval::operator()(Null_Ptr n)
  {
    return n;
  }

  Expression_Ptr Eval::operator()(Argument_Ptr a)
  {
    Expression_Obj val = a->value()->perform(this);
    bool is_rest_argument = a->is_rest_argument();
    bool is_keyword_argument = a->is_keyword_argument();

    if (a->is_rest_argument()) {
      if (val->concrete_type() == Expression::MAP) {
        is_rest_argument = false;
        is_keyword_argument = true;
      }
      else if(val->concrete_type() != Expression::LIST) {
        List_Obj wrapper = SASS_MEMORY_NEW(List,
                                        val->pstate(),
                                        0,
                                        SASS_COMMA,
                                        true);
        wrapper->append(val);
        val = wrapper;
      }
    }
    return SASS_MEMORY_NEW(Argument,
                           a->pstate(),
                           val,
                           a->name(),
                           is_rest_argument,
                           is_keyword_argument);
  }

  Expression_Ptr Eval::operator()(Arguments_Ptr a)
  {
    Arguments_Obj aa = SASS_MEMORY_NEW(Arguments, a->pstate());
    if (a->length() == 0) return aa.detach();
    for (size_t i = 0, L = a->length(); i < L; ++i) {
      Expression_Obj rv = (*a)[i]->perform(this);
      Argument_Ptr arg = Cast<Argument>(rv);
      if (!(arg->is_rest_argument() || arg->is_keyword_argument())) {
        aa->append(arg);
      }
    }

    if (a->has_rest_argument()) {
      Expression_Obj rest = a->get_rest_argument()->perform(this);
      Expression_Obj splat = Cast<Argument>(rest)->value()->perform(this);

      Sass_Separator separator = SASS_COMMA;
      List_Ptr ls = Cast<List>(splat);
      Map_Ptr ms = Cast<Map>(splat);

      List_Obj arglist = SASS_MEMORY_NEW(List,
                                      splat->pstate(),
                                      0,
                                      ls ? ls->separator() : separator,
                                      true);

      if (ls && ls->is_arglist()) {
        arglist->concat(ls);
      } else if (ms) {
        aa->append(SASS_MEMORY_NEW(Argument, splat->pstate(), ms, "", false, true));
      } else if (ls) {
        arglist->concat(ls);
      } else {
        arglist->append(splat);
      }
      if (arglist->length()) {
        aa->append(SASS_MEMORY_NEW(Argument, splat->pstate(), arglist, "", true));
      }
    }

    if (a->has_keyword_argument()) {
      Expression_Obj rv = a->get_keyword_argument()->perform(this);
      Argument_Ptr rvarg = Cast<Argument>(rv);
      Expression_Obj kwarg = rvarg->value()->perform(this);

      aa->append(SASS_MEMORY_NEW(Argument, kwarg->pstate(), kwarg, "", false, true));
    }
    return aa.detach();
  }

  Expression_Ptr Eval::operator()(Comment_Ptr c)
  {
    return 0;
  }

  Selector_List_Ptr Eval::operator()(Selector_List_Ptr s)
  {
    SelectorStack rv;
    Selector_List_Obj sl = SASS_MEMORY_NEW(Selector_List, s->pstate());
    sl->is_optional(s->is_optional());
    sl->media_block(s->media_block());
    sl->is_optional(s->is_optional());
    for (size_t i = 0, iL = s->length(); i < iL; ++i) {
      rv.push_back(operator()((*s)[i]));
    }

    // we should actually permutate parent first
    // but here we have permutated the selector first
    size_t round = 0;
    while (round != std::string::npos) {
      bool abort = true;
      for (size_t i = 0, iL = rv.size(); i < iL; ++i) {
        if (rv[i]->length() > round) {
          sl->append((*rv[i])[round]);
          abort = false;
        }
      }
      if (abort) {
        round = std::string::npos;
      } else {
        ++ round;
      }

    }
    return sl.detach();
  }


  Selector_List_Ptr Eval::operator()(Complex_Selector_Ptr s)
  {
    bool implicit_parent = !exp.old_at_root_without_rule;
    if (is_in_selector_schema) exp.selector_stack.push_back({});
    Selector_List_Obj resolved = s->resolve_parent_refs(exp.selector_stack, traces, implicit_parent);
    if (is_in_selector_schema) exp.selector_stack.pop_back();
    for (size_t i = 0; i < resolved->length(); i++) {
      Complex_Selector_Ptr is = resolved->at(i)->mutable_first();
      while (is) {
        if (is->head()) {
          is->head(operator()(is->head()));
        }
        is = is->tail();
      }
    }
    return resolved.detach();
  }

  Compound_Selector_Ptr Eval::operator()(Compound_Selector_Ptr s)
  {
    for (size_t i = 0; i < s->length(); i++) {
      Simple_Selector_Ptr ss = s->at(i);
      // skip parents here (called via resolve_parent_refs)
      if (ss == NULL || Cast<Parent_Selector>(ss)) continue;
      s->at(i) = Cast<Simple_Selector>(ss->perform(this));
    }
    return s;
  }

  Selector_List_Ptr Eval::operator()(Selector_Schema_Ptr s)
  {
    LOCAL_FLAG(is_in_selector_schema, true);
    // the parser will look for a brace to end the selector
    Expression_Obj sel = s->contents()->perform(this);
    std::string result_str(sel->to_string(options()));
    result_str = unquote(Util::rtrim(result_str));
    char* temp_cstr = sass_copy_c_string(result_str.c_str());
    ctx.strings.push_back(temp_cstr); // attach to context
    Parser p = Parser::from_c_str(temp_cstr, ctx, traces, s->pstate());
    p.last_media_block = s->media_block();
    // a selector schema may or may not connect to parent?
    bool chroot = s->connect_parent() == false;
    Selector_List_Obj sl = p.parse_selector_list(chroot);
    auto vec_str_rend = ctx.strings.rend();
    auto vec_str_rbegin = ctx.strings.rbegin();
    // remove the first item searching from the back
    // we cannot assume our item is still the last one
    // order is not important, so we can optimize this
    auto it = std::find(vec_str_rbegin, vec_str_rend, temp_cstr);
    // undefined behavior if not found!
    if (it != vec_str_rend) {
      // overwrite with last item
      *it = ctx.strings.back();
      // remove last one from vector
      ctx.strings.pop_back();
      // free temporary copy
      free(temp_cstr);
    }
    flag_is_in_selector_schema.reset();
    return operator()(sl);
  }

  Expression_Ptr Eval::operator()(Parent_Selector_Ptr p)
  {
    if (Selector_List_Obj pr = selector()) {
      exp.selector_stack.pop_back();
      Selector_List_Obj rv = operator()(pr);
      exp.selector_stack.push_back(rv);
      return rv.detach();
    } else {
      return SASS_MEMORY_NEW(Null, p->pstate());
    }
  }

  Expression_Ptr Eval::operator()(Parent_Reference_Ptr p)
  {
    if (Selector_List_Obj pr = selector()) {
      exp.selector_stack.pop_back();
      Selector_List_Obj rv = operator()(pr);
      exp.selector_stack.push_back(rv);
      return rv.detach();
    } else {
      return SASS_MEMORY_NEW(Null, p->pstate());
    }
  }

  Simple_Selector_Ptr Eval::operator()(Simple_Selector_Ptr s)
  {
    return s;
  }

  // hotfix to avoid invalid nested `:not` selectors
  // probably the wrong place, but this should ultimately
  // be fixed by implement superselector correctly for `:not`
  // first use of "find" (ATM only implemented for selectors)
  bool hasNotSelector(AST_Node_Obj obj) {
    if (Wrapped_Selector_Ptr w = Cast<Wrapped_Selector>(obj)) {
      return w->name() == ":not";
    }
    return false;
  }

  Wrapped_Selector_Ptr Eval::operator()(Wrapped_Selector_Ptr s)
  {

    if (s->name() == ":not") {
      if (exp.selector_stack.back()) {
        if (s->selector()->find(hasNotSelector)) {
          s->selector()->clear();
          s->name(" ");
        } else {
          for (size_t i = 0; i < s->selector()->length(); ++i) {
            Complex_Selector_Ptr cs = s->selector()->at(i);
            if (cs->tail()) {
              s->selector()->clear();
              s->name(" ");
            }
          }
        }
      }
    }
    return s;
  };

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.

#include <cctype>

namespace Sass {

  namespace Functions {

    bool string_argument(AST_Node_Obj obj) {
      String_Constant_Ptr s = Cast<String_Constant>(obj);
      if (s == nullptr) return false;
      const std::string& str = s->value();
      return starts_with(str, "calc(") ||
             starts_with(str, "var(");
    }

    void hsla_alpha_percent_deprecation(const ParserState& pstate, const std::string val)
    {

      std::string msg("Passing a percentage as the alpha value to hsla() will be interpreted");
      std::string tail("differently in future versions of Sass. For now, use " + val + " instead.");

      deprecated(msg, tail, false, pstate);

    }

    Signature rgb_sig = "rgb($red, $green, $blue)";
    BUILT_IN(rgb)
    {
      if (
        string_argument(env["$red"]) ||
        string_argument(env["$green"]) ||
        string_argument(env["$blue"])
      ) {
        return SASS_MEMORY_NEW(String_Constant, pstate, "rgb("
                                                        + env["$red"]->to_string()
                                                        + ", "
                                                        + env["$green"]->to_string()
                                                        + ", "
                                                        + env["$blue"]->to_string()
                                                        + ")"
        );
      }

      return SASS_MEMORY_NEW(Color_RGBA,
                             pstate,
                             COLOR_NUM("$red"),
                             COLOR_NUM("$green"),
                             COLOR_NUM("$blue"));
    }

    Signature rgba_4_sig = "rgba($red, $green, $blue, $alpha)";
    BUILT_IN(rgba_4)
    {
      if (
        string_argument(env["$red"]) ||
        string_argument(env["$green"]) ||
        string_argument(env["$blue"]) ||
        string_argument(env["$alpha"])
      ) {
        return SASS_MEMORY_NEW(String_Constant, pstate, "rgba("
                                                        + env["$red"]->to_string()
                                                        + ", "
                                                        + env["$green"]->to_string()
                                                        + ", "
                                                        + env["$blue"]->to_string()
                                                        + ", "
                                                        + env["$alpha"]->to_string()
                                                        + ")"
        );
      }

      return SASS_MEMORY_NEW(Color_RGBA,
                             pstate,
                             COLOR_NUM("$red"),
                             COLOR_NUM("$green"),
                             COLOR_NUM("$blue"),
                             ALPHA_NUM("$alpha"));
    }

    Signature rgba_2_sig = "rgba($color, $alpha)";
    BUILT_IN(rgba_2)
    {
      if (
        string_argument(env["$color"])
      ) {
        return SASS_MEMORY_NEW(String_Constant, pstate, "rgba("
                                                        + env["$color"]->to_string()
                                                        + ", "
                                                        + env["$alpha"]->to_string()
                                                        + ")"
        );
      }

      Color_RGBA_Obj c_arg = ARG("$color", Color)->toRGBA();

      if (
        string_argument(env["$alpha"])
      ) {
        std::stringstream strm;
        strm << "rgba("
                 << (int)c_arg->r() << ", "
                 << (int)c_arg->g() << ", "
                 << (int)c_arg->b() << ", "
                 << env["$alpha"]->to_string()
             << ")";
        return SASS_MEMORY_NEW(String_Constant, pstate, strm.str());
      }

      Color_RGBA_Obj new_c = SASS_MEMORY_COPY(c_arg);
      new_c->a(ALPHA_NUM("$alpha"));
      new_c->disp("");
      return new_c.detach();
    }

    ////////////////
    // RGB FUNCTIONS
    ////////////////

    Signature red_sig = "red($color)";
    BUILT_IN(red)
    {
      Color_RGBA_Obj color = ARG("$color", Color)->toRGBA();
      return SASS_MEMORY_NEW(Number, pstate, color->r());
    }

    Signature green_sig = "green($color)";
    BUILT_IN(green)
    {
      Color_RGBA_Obj color = ARG("$color", Color)->toRGBA();
      return SASS_MEMORY_NEW(Number, pstate, color->g());
    }

    Signature blue_sig = "blue($color)";
    BUILT_IN(blue)
    {
      Color_RGBA_Obj color = ARG("$color", Color)->toRGBA();
      return SASS_MEMORY_NEW(Number, pstate, color->b());
    }

    Color_RGBA* colormix(Context& ctx, ParserState& pstate, Color* color1, Color* color2, double weight) {
      Color_RGBA_Obj c1 = color1->toRGBA();
      Color_RGBA_Obj c2 = color2->toRGBA();
      double p = weight/100;
      double w = 2*p - 1;
      double a = c1->a() - c2->a();

      double w1 = (((w * a == -1) ? w : (w + a)/(1 + w*a)) + 1)/2.0;
      double w2 = 1 - w1;

      return SASS_MEMORY_NEW(Color_RGBA,
                             pstate,
                             Sass::round(w1*c1->r() + w2*c2->r(), ctx.c_options.precision),
                             Sass::round(w1*c1->g() + w2*c2->g(), ctx.c_options.precision),
                             Sass::round(w1*c1->b() + w2*c2->b(), ctx.c_options.precision),
                             c1->a()*p + c2->a()*(1-p));
    }

    Signature mix_sig = "mix($color-1, $color-2, $weight: 50%)";
    BUILT_IN(mix)
    {
      Color_Obj  color1 = ARG("$color-1", Color);
      Color_Obj  color2 = ARG("$color-2", Color);
      double weight = DARG_U_PRCT("$weight");
      return colormix(ctx, pstate, color1, color2, weight);

    }

    ////////////////
    // HSL FUNCTIONS
    ////////////////

    Signature hsl_sig = "hsl($hue, $saturation, $lightness)";
    BUILT_IN(hsl)
    {
      if (
        string_argument(env["$hue"]) ||
        string_argument(env["$saturation"]) ||
        string_argument(env["$lightness"])
      ) {
        return SASS_MEMORY_NEW(String_Constant, pstate, "hsl("
                                                        + env["$hue"]->to_string()
                                                        + ", "
                                                        + env["$saturation"]->to_string()
                                                        + ", "
                                                        + env["$lightness"]->to_string()
                                                        + ")"
        );
      }

      return SASS_MEMORY_NEW(Color_HSLA,
        pstate,
        ARGVAL("$hue"),
        ARGVAL("$saturation"),
        ARGVAL("$lightness"),
        1.0);

    }

    Signature hsla_sig = "hsla($hue, $saturation, $lightness, $alpha)";
    BUILT_IN(hsla)
    {
      if (
        string_argument(env["$hue"]) ||
        string_argument(env["$saturation"]) ||
        string_argument(env["$lightness"]) ||
        string_argument(env["$alpha"])
      ) {
        return SASS_MEMORY_NEW(String_Constant, pstate, "hsla("
                                                        + env["$hue"]->to_string()
                                                        + ", "
                                                        + env["$saturation"]->to_string()
                                                        + ", "
                                                        + env["$lightness"]->to_string()
                                                        + ", "
                                                        + env["$alpha"]->to_string()
                                                        + ")"
        );
      }

      Number_Ptr alpha = ARG("$alpha", Number);
      if (alpha && alpha->unit() == "%") {
        Number_Obj val = SASS_MEMORY_COPY(alpha);
        val->numerators.clear(); // convert
        val->value(val->value() / 100.0);
        std::string nr(val->to_string(ctx.c_options));
        hsla_alpha_percent_deprecation(pstate, nr);
      }

      return SASS_MEMORY_NEW(Color_HSLA,
        pstate,
        ARGVAL("$hue"),
        ARGVAL("$saturation"),
        ARGVAL("$lightness"),
        ARGVAL("$alpha"));

    }

    /////////////////////////////////////////////////////////////////////////
    // Query functions
    /////////////////////////////////////////////////////////////////////////

    Signature hue_sig = "hue($color)";
    BUILT_IN(hue)
    {
      Color_HSLA_Obj col = ARG("$color", Color)->toHSLA();
      return SASS_MEMORY_NEW(Number, pstate, col->h(), "deg");
    }

    Signature saturation_sig = "saturation($color)";
    BUILT_IN(saturation)
    {
      Color_HSLA_Obj col = ARG("$color", Color)->toHSLA();
      return SASS_MEMORY_NEW(Number, pstate, col->s(), "%");
    }

    Signature lightness_sig = "lightness($color)";
    BUILT_IN(lightness)
    {
      Color_HSLA_Obj col = ARG("$color", Color)->toHSLA();
      return SASS_MEMORY_NEW(Number, pstate, col->l(), "%");
    }

    /////////////////////////////////////////////////////////////////////////
    // HSL manipulation functions
    /////////////////////////////////////////////////////////////////////////

    Signature adjust_hue_sig = "adjust-hue($color, $degrees)";
    BUILT_IN(adjust_hue)
    {
      Color_Ptr col = ARG("$color", Color);
      double degrees = ARGVAL("$degrees");
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->h(absmod(copy->h() + degrees, 360.0));
      return copy.detach();
    }

    Signature lighten_sig = "lighten($color, $amount)";
    BUILT_IN(lighten)
    {
      Color_Ptr col = ARG("$color", Color);
      double amount = DARG_U_PRCT("$amount");
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->l(clip(copy->l() + amount, 0.0, 100.0));
      return copy.detach();

    }

    Signature darken_sig = "darken($color, $amount)";
    BUILT_IN(darken)
    {
      Color_Ptr col = ARG("$color", Color);
      double amount = DARG_U_PRCT("$amount");
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->l(clip(copy->l() - amount, 0.0, 100.0));
      return copy.detach();
    }

    Signature saturate_sig = "saturate($color, $amount: false)";
    BUILT_IN(saturate)
    {
      // CSS3 filter function overload: pass literal through directly
      if (!Cast<Number>(env["$amount"])) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "saturate(" + env["$color"]->to_string(ctx.c_options) + ")");
      }

      Color_Ptr col = ARG("$color", Color);
      double amount = DARG_U_PRCT("$amount");
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->s(clip(copy->s() + amount, 0.0, 100.0));
      return copy.detach();
    }

    Signature desaturate_sig = "desaturate($color, $amount)";
    BUILT_IN(desaturate)
    {
      Color_Ptr col = ARG("$color", Color);
      double amount = DARG_U_PRCT("$amount");
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->s(clip(copy->s() - amount, 0.0, 100.0));
      return copy.detach();
    }

    Signature grayscale_sig = "grayscale($color)";
    BUILT_IN(grayscale)
    {
      // CSS3 filter function overload: pass literal through directly
      Number_Ptr amount = Cast<Number>(env["$color"]);
      if (amount) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "grayscale(" + amount->to_string(ctx.c_options) + ")");
      }

      Color_Ptr col = ARG("$color", Color);
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->s(0.0); // just reset saturation
      return copy.detach();
    }

    /////////////////////////////////////////////////////////////////////////
    // Misc manipulation functions
    /////////////////////////////////////////////////////////////////////////

    Signature complement_sig = "complement($color)";
    BUILT_IN(complement)
    {
      Color_Ptr col = ARG("$color", Color);
      Color_HSLA_Obj copy = col->toHSLA(true);
      copy->h(absmod(copy->h() - 180.0, 360.0));
      return copy.detach();
    }

    Signature invert_sig = "invert($color, $weight: 100%)";
    BUILT_IN(invert)
    {
      // CSS3 filter function overload: pass literal through directly
      Number_Ptr amount = Cast<Number>(env["$color"]);
      if (amount) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "invert(" + amount->to_string(ctx.c_options) + ")");
      }

      Color_Ptr col = ARG("$color", Color);
      double weight = DARG_U_PRCT("$weight");
      Color_RGBA_Obj inv = col->toRGBA(true);
      inv->r(clip(255.0 - inv->r(), 0.0, 255.0));
      inv->g(clip(255.0 - inv->g(), 0.0, 255.0));
      inv->b(clip(255.0 - inv->b(), 0.0, 255.0));
      return colormix(ctx, pstate, inv, col, weight);
    }

    /////////////////////////////////////////////////////////////////////////
    // Opacity functions
    /////////////////////////////////////////////////////////////////////////

    Signature alpha_sig = "alpha($color)";
    Signature opacity_sig = "opacity($color)";
    BUILT_IN(alpha)
    {
      String_Constant_Ptr ie_kwd = Cast<String_Constant>(env["$color"]);
      if (ie_kwd) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "alpha(" + ie_kwd->value() + ")");
      }

      // CSS3 filter function overload: pass literal through directly
      Number_Ptr amount = Cast<Number>(env["$color"]);
      if (amount) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "opacity(" + amount->to_string(ctx.c_options) + ")");
      }

      return SASS_MEMORY_NEW(Number, pstate, ARG("$color", Color)->a());
    }

    Signature opacify_sig = "opacify($color, $amount)";
    Signature fade_in_sig = "fade-in($color, $amount)";
    BUILT_IN(opacify)
    {
      Color_Ptr col = ARG("$color", Color);
      double amount = DARG_U_FACT("$amount");
      Color_Obj copy = SASS_MEMORY_COPY(col);
      copy->a(clip(col->a() + amount, 0.0, 1.0));
      return copy.detach();
    }

    Signature transparentize_sig = "transparentize($color, $amount)";
    Signature fade_out_sig = "fade-out($color, $amount)";
    BUILT_IN(transparentize)
    {
      Color_Ptr col = ARG("$color", Color);
      double amount = DARG_U_FACT("$amount");
      Color_Obj copy = SASS_MEMORY_COPY(col);
      copy->a(std::max(col->a() - amount, 0.0));
      return copy.detach();
    }

    ////////////////////////
    // OTHER COLOR FUNCTIONS
    ////////////////////////

    Signature adjust_color_sig = "adjust-color($color, $red: false, $green: false, $blue: false, $hue: false, $saturation: false, $lightness: false, $alpha: false)";
    BUILT_IN(adjust_color)
    {
      Color_Ptr col = ARG("$color", Color);
      Number_Ptr r = Cast<Number>(env["$red"]);
      Number_Ptr g = Cast<Number>(env["$green"]);
      Number_Ptr b = Cast<Number>(env["$blue"]);
      Number_Ptr h = Cast<Number>(env["$hue"]);
      Number_Ptr s = Cast<Number>(env["$saturation"]);
      Number_Ptr l = Cast<Number>(env["$lightness"]);
      Number_Ptr a = Cast<Number>(env["$alpha"]);

      bool rgb = r || g || b;
      bool hsl = h || s || l;

      if (rgb && hsl) {
        error("Cannot specify HSL and RGB values for a color at the same time for `adjust-color'", pstate, traces);
      }
      else if (rgb) {
        Color_RGBA_Obj c = col->toRGBA(true);
        if (r) c->r(c->r() + DARG_R_BYTE("$red"));
        if (g) c->g(c->g() + DARG_R_BYTE("$green"));
        if (b) c->b(c->b() + DARG_R_BYTE("$blue"));
        if (a) c->a(c->a() + DARG_R_FACT("$alpha"));
        return c.detach();
      }
      else if (hsl) {
        Color_HSLA_Obj c = col->toHSLA(true);
        if (h) c->h(c->h() + absmod(h->value(), 360.0));
        if (s) c->s(c->s() + DARG_R_PRCT("$saturation"));
        if (l) c->l(c->l() + DARG_R_PRCT("$lightness"));
        if (a) c->a(c->a() + DARG_R_FACT("$alpha"));
        return c.detach();
      }
      else if (a) {
        Color_Obj c = SASS_MEMORY_COPY(col);
        c->a(c->a() + DARG_R_FACT("$alpha"));
        c->a(clip(c->a(), 0.0, 1.0));
        return c.detach();
      }
      error("not enough arguments for `adjust-color'", pstate, traces);
      // unreachable
      return col;
    }

    Signature scale_color_sig = "scale-color($color, $red: false, $green: false, $blue: false, $hue: false, $saturation: false, $lightness: false, $alpha: false)";
    BUILT_IN(scale_color)
    {
      Color_Ptr col = ARG("$color", Color);
      Number_Ptr r = Cast<Number>(env["$red"]);
      Number_Ptr g = Cast<Number>(env["$green"]);
      Number_Ptr b = Cast<Number>(env["$blue"]);
      Number_Ptr h = Cast<Number>(env["$hue"]);
      Number_Ptr s = Cast<Number>(env["$saturation"]);
      Number_Ptr l = Cast<Number>(env["$lightness"]);
      Number_Ptr a = Cast<Number>(env["$alpha"]);

      bool rgb = r || g || b;
      bool hsl = h || s || l;

      if (rgb && hsl) {
        error("Cannot specify HSL and RGB values for a color at the same time for `scale-color'", pstate, traces);
      }
      else if (rgb) {
        Color_RGBA_Obj c = col->toRGBA(true);
        double rscale = (r ? DARG_R_PRCT("$red") : 0.0) / 100.0;
        double gscale = (g ? DARG_R_PRCT("$green") : 0.0) / 100.0;
        double bscale = (b ? DARG_R_PRCT("$blue") : 0.0) / 100.0;
        double ascale = (a ? DARG_R_PRCT("$alpha") : 0.0) / 100.0;
        if (rscale) c->r(c->r() + rscale * (rscale > 0.0 ? 255.0 - c->r() : c->r()));
        if (gscale) c->g(c->g() + gscale * (gscale > 0.0 ? 255.0 - c->g() : c->g()));
        if (bscale) c->b(c->b() + bscale * (bscale > 0.0 ? 255.0 - c->b() : c->b()));
        if (ascale) c->a(c->a() + ascale * (ascale > 0.0 ? 1.0 - c->a() : c->a()));
        return c.detach();
      }
      else if (hsl) {
        Color_HSLA_Obj c = col->toHSLA(true);
        double hscale = (h ? DARG_R_PRCT("$hue") : 0.0) / 100.0;
        double sscale = (s ? DARG_R_PRCT("$saturation") : 0.0) / 100.0;
        double lscale = (l ? DARG_R_PRCT("$lightness") : 0.0) / 100.0;
        double ascale = (a ? DARG_R_PRCT("$alpha") : 0.0) / 100.0;
        if (hscale) c->h(c->h() + hscale * (hscale > 0.0 ? 360.0 - c->h() : c->h()));
        if (sscale) c->s(c->s() + sscale * (sscale > 0.0 ? 100.0 - c->l() : c->s()));
        if (lscale) c->l(c->l() + lscale * (lscale > 0.0 ? 100.0 - c->l() : c->l()));
        if (ascale) c->a(c->a() + ascale * (ascale > 0.0 ? 1.0 - c->a() : c->a()));
        return c.detach();
      }
      else if (a) {
        Color_Obj c = SASS_MEMORY_COPY(col);
        double ascale = DARG_R_PRCT("$alpha") / 100.0;
        c->a(c->a() + ascale * (ascale > 0.0 ? 1.0 - c->a() : c->a()));
        c->a(clip(c->a(), 0.0, 1.0));
        return c.detach();
      }
      error("not enough arguments for `scale-color'", pstate, traces);
      // unreachable
      return col;
    }

    Signature change_color_sig = "change-color($color, $red: false, $green: false, $blue: false, $hue: false, $saturation: false, $lightness: false, $alpha: false)";
    BUILT_IN(change_color)
    {
      Color_Ptr col = ARG("$color", Color);
      Number_Ptr r = Cast<Number>(env["$red"]);
      Number_Ptr g = Cast<Number>(env["$green"]);
      Number_Ptr b = Cast<Number>(env["$blue"]);
      Number_Ptr h = Cast<Number>(env["$hue"]);
      Number_Ptr s = Cast<Number>(env["$saturation"]);
      Number_Ptr l = Cast<Number>(env["$lightness"]);
      Number_Ptr a = Cast<Number>(env["$alpha"]);

      bool rgb = r || g || b;
      bool hsl = h || s || l;

      if (rgb && hsl) {
        error("Cannot specify HSL and RGB values for a color at the same time for `change-color'", pstate, traces);
      }
      else if (rgb) {
        Color_RGBA_Obj c = col->toRGBA(true);
        if (r) c->r(DARG_U_BYTE("$red"));
        if (g) c->g(DARG_U_BYTE("$green"));
        if (b) c->b(DARG_U_BYTE("$blue"));
        if (a) c->a(DARG_U_FACT("$alpha"));
        return c.detach();
      }
      else if (hsl) {
        Color_HSLA_Obj c = col->toHSLA(true);
        if (h) c->h(absmod(h->value(), 360.0));
        if (s) c->s(DARG_U_PRCT("$saturation"));
        if (l) c->l(DARG_U_PRCT("$lightness"));
        if (a) c->a(DARG_U_FACT("$alpha"));
        return c.detach();
      }
      else if (a) {
        Color_Obj c = SASS_MEMORY_COPY(col);
        c->a(clip(DARG_U_FACT("$alpha"), 0.0, 1.0));
        return c.detach();
      }
      error("not enough arguments for `change-color'", pstate, traces);
      // unreachable
      return col;
    }

    Signature ie_hex_str_sig = "ie-hex-str($color)";
    BUILT_IN(ie_hex_str)
    {
      Color_Ptr col = ARG("$color", Color);
      Color_RGBA_Obj c = col->toRGBA();
      double r = clip(c->r(), 0.0, 255.0);
      double g = clip(c->g(), 0.0, 255.0);
      double b = clip(c->b(), 0.0, 255.0);
      double a = clip(c->a(), 0.0, 1.0) * 255.0;

      std::stringstream ss;
      ss << '#' << std::setw(2) << std::setfill('0');
      ss << std::hex << std::setw(2) << static_cast<unsigned long>(Sass::round(a, ctx.c_options.precision));
      ss << std::hex << std::setw(2) << static_cast<unsigned long>(Sass::round(r, ctx.c_options.precision));
      ss << std::hex << std::setw(2) << static_cast<unsigned long>(Sass::round(g, ctx.c_options.precision));
      ss << std::hex << std::setw(2) << static_cast<unsigned long>(Sass::round(b, ctx.c_options.precision));

      std::string result(ss.str());
      for (size_t i = 0, L = result.length(); i < L; ++i) {
        result[i] = std::toupper(result[i]);
      }
      return SASS_MEMORY_NEW(String_Quoted, pstate, result);
    }

  }

}

namespace Sass {

  void bind(std::string type, std::string name, Parameters_Obj ps, Arguments_Obj as, Env* env, Eval* eval, Backtraces& traces)
  {
    std::string callee(type + " " + name);

    std::map<std::string, Parameter_Obj> param_map;
    List_Obj varargs = SASS_MEMORY_NEW(List, as->pstate());
    varargs->is_arglist(true); // enable keyword size handling

    for (size_t i = 0, L = as->length(); i < L; ++i) {
      if (auto str = Cast<String_Quoted>((*as)[i]->value())) {
        // force optional quotes (only if needed)
        if (str->quote_mark()) {
          str->quote_mark('*');
        }
      }
    }

    // Set up a map to ensure named arguments refer to actual parameters. Also
    // eval each default value left-to-right, wrt env, populating env as we go.
    for (size_t i = 0, L = ps->length(); i < L; ++i) {
      Parameter_Obj  p = ps->at(i);
      param_map[p->name()] = p;
      // if (p->default_value()) {
      //   env->local_frame()[p->name()] = p->default_value()->perform(eval->with(env));
      // }
    }

    // plug in all args; if we have leftover params, deal with it later
    size_t ip = 0, LP = ps->length();
    size_t ia = 0, LA = as->length();
    while (ia < LA) {
      Argument_Obj a = as->at(ia);
      if (ip >= LP) {
        // skip empty rest arguments
        if (a->is_rest_argument()) {
          if (List_Obj l = Cast<List>(a->value())) {
            if (l->length() == 0) {
              ++ ia; continue;
            }
          }
        }
        std::stringstream msg;
        msg << "wrong number of arguments (" << LA << " for " << LP << ")";
        msg << " for `" << name << "'";
        return error(msg.str(), as->pstate(), traces);
      }
      Parameter_Obj p = ps->at(ip);

      // If the current parameter is the rest parameter, process and break the loop
      if (p->is_rest_parameter()) {
        // The next argument by coincidence provides a rest argument
        if (a->is_rest_argument()) {

          // We should always get a list for rest arguments
          if (List_Obj rest = Cast<List>(a->value())) {
              // create a new list object for wrapped items
              List_Ptr arglist = SASS_MEMORY_NEW(List,
                                              p->pstate(),
                                              0,
                                              rest->separator(),
                                              true);
              // wrap each item from list as an argument
              for (Expression_Obj item : rest->elements()) {
                if (Argument_Obj arg = Cast<Argument>(item)) {
                  arglist->append(SASS_MEMORY_COPY(arg)); // copy
                } else {
                  arglist->append(SASS_MEMORY_NEW(Argument,
                                                  item->pstate(),
                                                  item,
                                                  "",
                                                  false,
                                                  false));
                }
              }
              // assign new arglist to environment
              env->local_frame()[p->name()] = arglist;
            }
          // invalid state
          else {
            throw std::runtime_error("invalid state");
          }
        } else if (a->is_keyword_argument()) {

          // expand keyword arguments into their parameters
          List_Ptr arglist = SASS_MEMORY_NEW(List, p->pstate(), 0, SASS_COMMA, true);
          env->local_frame()[p->name()] = arglist;
          Map_Obj argmap = Cast<Map>(a->value());
          for (auto key : argmap->keys()) {
            if (String_Constant_Obj str = Cast<String_Constant>(key)) {
              std::string param = unquote(str->value());
              arglist->append(SASS_MEMORY_NEW(Argument,
                                              key->pstate(),
                                              argmap->at(key),
                                              "$" + param,
                                              false,
                                              false));
            } else {
              traces.push_back(Backtrace(key->pstate()));
              throw Exception::InvalidVarKwdType(key->pstate(), traces, key->inspect(), a);
            }
          }

        } else {

          // create a new list object for wrapped items
          List_Obj arglist = SASS_MEMORY_NEW(List,
                                          p->pstate(),
                                          0,
                                          SASS_COMMA,
                                          true);
          // consume the next args
          while (ia < LA) {
            // get and post inc
            a = (*as)[ia++];
            // maybe we have another list as argument
            List_Obj ls = Cast<List>(a->value());
            // skip any list completely if empty
            if (ls && ls->empty() && a->is_rest_argument()) continue;

            Expression_Obj value = a->value();
            if (Argument_Obj arg = Cast<Argument>(value)) {
              arglist->append(arg);
            }
            // check if we have rest argument
            else if (a->is_rest_argument()) {
              // preserve the list separator from rest args
              if (List_Obj rest = Cast<List>(a->value())) {
                arglist->separator(rest->separator());

                for (size_t i = 0, L = rest->length(); i < L; ++i) {
                  Expression_Obj obj = rest->value_at_index(i);
                  arglist->append(SASS_MEMORY_NEW(Argument,
                                                obj->pstate(),
                                                obj,
                                                "",
                                                false,
                                                false));
                }
              }
              // no more arguments
              break;
            }
            // wrap all other value types into Argument
            else {
              arglist->append(SASS_MEMORY_NEW(Argument,
                                            a->pstate(),
                                            a->value(),
                                            a->name(),
                                            false,
                                            false));
            }
          }
          // assign new arglist to environment
          env->local_frame()[p->name()] = arglist;
        }
        // consumed parameter
        ++ip;
        // no more paramaters
        break;
      }

      // If the current argument is the rest argument, extract a value for processing
      else if (a->is_rest_argument()) {
        // normal param and rest arg
        List_Obj arglist = Cast<List>(a->value());
        if (!arglist) {
          if (Expression_Obj arg = Cast<Expression>(a->value())) {
            arglist = SASS_MEMORY_NEW(List, a->pstate(), 1);
            arglist->append(arg);
          }
        }

        // empty rest arg - treat all args as default values
        if (!arglist || !arglist->length()) {
          break;
        } else {
          if (arglist->length() > LP - ip && !ps->has_rest_parameter()) {
            size_t arg_count = (arglist->length() + LA - 1);
            std::stringstream msg;
            msg << callee << " takes " << LP;
            msg << (LP == 1 ? " argument" : " arguments");
            msg << " but " << arg_count;
            msg << (arg_count == 1 ? " was passed" : " were passed.");
            deprecated_bind(msg.str(), as->pstate());

            while (arglist->length() > LP - ip) {
              arglist->elements().erase(arglist->elements().end() - 1);
            }
          }
        }
        // otherwise move one of the rest args into the param, converting to argument if necessary
        Expression_Obj obj = arglist->at(0);
        if (!(a = Cast<Argument>(obj))) {
          Expression_Ptr a_to_convert = obj;
          a = SASS_MEMORY_NEW(Argument,
                              a_to_convert->pstate(),
                              a_to_convert,
                              "",
                              false,
                              false);
        }
        arglist->elements().erase(arglist->elements().begin());
        if (!arglist->length() || (!arglist->is_arglist() && ip + 1 == LP)) {
          ++ia;
        }

      } else if (a->is_keyword_argument()) {
        Map_Obj argmap = Cast<Map>(a->value());

        for (auto key : argmap->keys()) {
          String_Constant_Ptr val = Cast<String_Constant>(key);
          if (val == NULL) {
            traces.push_back(Backtrace(key->pstate()));
            throw Exception::InvalidVarKwdType(key->pstate(), traces, key->inspect(), a);
          }
          std::string param = "$" + unquote(val->value());

          if (!param_map.count(param)) {
            std::stringstream msg;
            msg << callee << " has no parameter named " << param;
            error(msg.str(), a->pstate(), traces);
          }
          env->local_frame()[param] = argmap->at(key);
        }
        ++ia;
        continue;
      } else {
        ++ia;
      }

      if (a->name().empty()) {
        if (env->has_local(p->name())) {
          std::stringstream msg;
          msg << "parameter " << p->name()
          << " provided more than once in call to " << callee;
          error(msg.str(), a->pstate(), traces);
        }
        // ordinal arg -- bind it to the next param
        env->local_frame()[p->name()] = a->value();
        ++ip;
      }
      else {
        // named arg -- bind it to the appropriately named param
        if (!param_map.count(a->name())) {
          if (ps->has_rest_parameter()) {
            varargs->append(a);
          } else {
            std::stringstream msg;
            msg << callee << " has no parameter named " << a->name();
            error(msg.str(), a->pstate(), traces);
          }
        }
        if (param_map[a->name()]) {
          if (param_map[a->name()]->is_rest_parameter()) {
            std::stringstream msg;
            msg << "argument " << a->name() << " of " << callee
                << "cannot be used as named argument";
            error(msg.str(), a->pstate(), traces);
          }
        }
        if (env->has_local(a->name())) {
          std::stringstream msg;
          msg << "parameter " << p->name()
              << "provided more than once in call to " << callee;
          error(msg.str(), a->pstate(), traces);
        }
        env->local_frame()[a->name()] = a->value();
      }
    }
    // EO while ia

    // If we make it here, we're out of args but may have leftover params.
    // That's only okay if they have default values, or were already bound by
    // named arguments, or if it's a single rest-param.
    for (size_t i = ip; i < LP; ++i) {
      Parameter_Obj leftover = ps->at(i);
      // cerr << "env for default params:" << endl;
      // env->print();
      // cerr << "********" << endl;
      if (!env->has_local(leftover->name())) {
        if (leftover->is_rest_parameter()) {
          env->local_frame()[leftover->name()] = varargs;
        }
        else if (leftover->default_value()) {
          Expression_Ptr dv = leftover->default_value()->perform(eval);
          env->local_frame()[leftover->name()] = dv;
        }
        else {
          // param is unbound and has no default value -- error
          throw Exception::MissingArgument(as->pstate(), traces, name, leftover->name(), type);
        }
      }
    }

    return;
  }


}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  namespace ColorNames
  {
    const char aliceblue [] = "aliceblue";
    const char antiquewhite [] = "antiquewhite";
    const char cyan [] = "cyan";
    const char aqua [] = "aqua";
    const char aquamarine [] = "aquamarine";
    const char azure [] = "azure";
    const char beige [] = "beige";
    const char bisque [] = "bisque";
    const char black [] = "black";
    const char blanchedalmond [] = "blanchedalmond";
    const char blue [] = "blue";
    const char blueviolet [] = "blueviolet";
    const char brown [] = "brown";
    const char burlywood [] = "burlywood";
    const char cadetblue [] = "cadetblue";
    const char chartreuse [] = "chartreuse";
    const char chocolate [] = "chocolate";
    const char coral [] = "coral";
    const char cornflowerblue [] = "cornflowerblue";
    const char cornsilk [] = "cornsilk";
    const char crimson [] = "crimson";
    const char darkblue [] = "darkblue";
    const char darkcyan [] = "darkcyan";
    const char darkgoldenrod [] = "darkgoldenrod";
    const char darkgray [] = "darkgray";
    const char darkgrey [] = "darkgrey";
    const char darkgreen [] = "darkgreen";
    const char darkkhaki [] = "darkkhaki";
    const char darkmagenta [] = "darkmagenta";
    const char darkolivegreen [] = "darkolivegreen";
    const char darkorange [] = "darkorange";
    const char darkorchid [] = "darkorchid";
    const char darkred [] = "darkred";
    const char darksalmon [] = "darksalmon";
    const char darkseagreen [] = "darkseagreen";
    const char darkslateblue [] = "darkslateblue";
    const char darkslategray [] = "darkslategray";
    const char darkslategrey [] = "darkslategrey";
    const char darkturquoise [] = "darkturquoise";
    const char darkviolet [] = "darkviolet";
    const char deeppink [] = "deeppink";
    const char deepskyblue [] = "deepskyblue";
    const char dimgray [] = "dimgray";
    const char dimgrey [] = "dimgrey";
    const char dodgerblue [] = "dodgerblue";
    const char firebrick [] = "firebrick";
    const char floralwhite [] = "floralwhite";
    const char forestgreen [] = "forestgreen";
    const char magenta [] = "magenta";
    const char fuchsia [] = "fuchsia";
    const char gainsboro [] = "gainsboro";
    const char ghostwhite [] = "ghostwhite";
    const char gold [] = "gold";
    const char goldenrod [] = "goldenrod";
    const char gray [] = "gray";
    const char grey [] = "grey";
    const char green [] = "green";
    const char greenyellow [] = "greenyellow";
    const char honeydew [] = "honeydew";
    const char hotpink [] = "hotpink";
    const char indianred [] = "indianred";
    const char indigo [] = "indigo";
    const char ivory [] = "ivory";
    const char khaki [] = "khaki";
    const char lavender [] = "lavender";
    const char lavenderblush [] = "lavenderblush";
    const char lawngreen [] = "lawngreen";
    const char lemonchiffon [] = "lemonchiffon";
    const char lightblue [] = "lightblue";
    const char lightcoral [] = "lightcoral";
    const char lightcyan [] = "lightcyan";
    const char lightgoldenrodyellow [] = "lightgoldenrodyellow";
    const char lightgray [] = "lightgray";
    const char lightgrey [] = "lightgrey";
    const char lightgreen [] = "lightgreen";
    const char lightpink [] = "lightpink";
    const char lightsalmon [] = "lightsalmon";
    const char lightseagreen [] = "lightseagreen";
    const char lightskyblue [] = "lightskyblue";
    const char lightslategray [] = "lightslategray";
    const char lightslategrey [] = "lightslategrey";
    const char lightsteelblue [] = "lightsteelblue";
    const char lightyellow [] = "lightyellow";
    const char lime [] = "lime";
    const char limegreen [] = "limegreen";
    const char linen [] = "linen";
    const char maroon [] = "maroon";
    const char mediumaquamarine [] = "mediumaquamarine";
    const char mediumblue [] = "mediumblue";
    const char mediumorchid [] = "mediumorchid";
    const char mediumpurple [] = "mediumpurple";
    const char mediumseagreen [] = "mediumseagreen";
    const char mediumslateblue [] = "mediumslateblue";
    const char mediumspringgreen [] = "mediumspringgreen";
    const char mediumturquoise [] = "mediumturquoise";
    const char mediumvioletred [] = "mediumvioletred";
    const char midnightblue [] = "midnightblue";
    const char mintcream [] = "mintcream";
    const char mistyrose [] = "mistyrose";
    const char moccasin [] = "moccasin";
    const char navajowhite [] = "navajowhite";
    const char navy [] = "navy";
    const char oldlace [] = "oldlace";
    const char olive [] = "olive";
    const char olivedrab [] = "olivedrab";
    const char orange [] = "orange";
    const char orangered [] = "orangered";
    const char orchid [] = "orchid";
    const char palegoldenrod [] = "palegoldenrod";
    const char palegreen [] = "palegreen";
    const char paleturquoise [] = "paleturquoise";
    const char palevioletred [] = "palevioletred";
    const char papayawhip [] = "papayawhip";
    const char peachpuff [] = "peachpuff";
    const char peru [] = "peru";
    const char pink [] = "pink";
    const char plum [] = "plum";
    const char powderblue [] = "powderblue";
    const char purple [] = "purple";
    const char red [] = "red";
    const char rosybrown [] = "rosybrown";
    const char royalblue [] = "royalblue";
    const char saddlebrown [] = "saddlebrown";
    const char salmon [] = "salmon";
    const char sandybrown [] = "sandybrown";
    const char seagreen [] = "seagreen";
    const char seashell [] = "seashell";
    const char sienna [] = "sienna";
    const char silver [] = "silver";
    const char skyblue [] = "skyblue";
    const char slateblue [] = "slateblue";
    const char slategray [] = "slategray";
    const char slategrey [] = "slategrey";
    const char snow [] = "snow";
    const char springgreen [] = "springgreen";
    const char steelblue [] = "steelblue";
    const char tan [] = "tan";
    const char teal [] = "teal";
    const char thistle [] = "thistle";
    const char tomato [] = "tomato";
    const char turquoise [] = "turquoise";
    const char violet [] = "violet";
    const char wheat [] = "wheat";
    const char white [] = "white";
    const char whitesmoke [] = "whitesmoke";
    const char yellow [] = "yellow";
    const char yellowgreen [] = "yellowgreen";
    const char rebeccapurple [] = "rebeccapurple";
    const char transparent [] = "transparent";
  }

  namespace Colors {
    const ParserState color_table("[COLOR TABLE]");
    const Color_RGBA aliceblue(color_table, 240, 248, 255, 1);
    const Color_RGBA antiquewhite(color_table, 250, 235, 215, 1);
    const Color_RGBA cyan(color_table, 0, 255, 255, 1);
    const Color_RGBA aqua(color_table, 0, 255, 255, 1);
    const Color_RGBA aquamarine(color_table, 127, 255, 212, 1);
    const Color_RGBA azure(color_table, 240, 255, 255, 1);
    const Color_RGBA beige(color_table, 245, 245, 220, 1);
    const Color_RGBA bisque(color_table, 255, 228, 196, 1);
    const Color_RGBA black(color_table, 0, 0, 0, 1);
    const Color_RGBA blanchedalmond(color_table, 255, 235, 205, 1);
    const Color_RGBA blue(color_table, 0, 0, 255, 1);
    const Color_RGBA blueviolet(color_table, 138, 43, 226, 1);
    const Color_RGBA brown(color_table, 165, 42, 42, 1);
    const Color_RGBA burlywood(color_table, 222, 184, 135, 1);
    const Color_RGBA cadetblue(color_table, 95, 158, 160, 1);
    const Color_RGBA chartreuse(color_table, 127, 255, 0, 1);
    const Color_RGBA chocolate(color_table, 210, 105, 30, 1);
    const Color_RGBA coral(color_table, 255, 127, 80, 1);
    const Color_RGBA cornflowerblue(color_table, 100, 149, 237, 1);
    const Color_RGBA cornsilk(color_table, 255, 248, 220, 1);
    const Color_RGBA crimson(color_table, 220, 20, 60, 1);
    const Color_RGBA darkblue(color_table, 0, 0, 139, 1);
    const Color_RGBA darkcyan(color_table, 0, 139, 139, 1);
    const Color_RGBA darkgoldenrod(color_table, 184, 134, 11, 1);
    const Color_RGBA darkgray(color_table, 169, 169, 169, 1);
    const Color_RGBA darkgrey(color_table, 169, 169, 169, 1);
    const Color_RGBA darkgreen(color_table, 0, 100, 0, 1);
    const Color_RGBA darkkhaki(color_table, 189, 183, 107, 1);
    const Color_RGBA darkmagenta(color_table, 139, 0, 139, 1);
    const Color_RGBA darkolivegreen(color_table, 85, 107, 47, 1);
    const Color_RGBA darkorange(color_table, 255, 140, 0, 1);
    const Color_RGBA darkorchid(color_table, 153, 50, 204, 1);
    const Color_RGBA darkred(color_table, 139, 0, 0, 1);
    const Color_RGBA darksalmon(color_table, 233, 150, 122, 1);
    const Color_RGBA darkseagreen(color_table, 143, 188, 143, 1);
    const Color_RGBA darkslateblue(color_table, 72, 61, 139, 1);
    const Color_RGBA darkslategray(color_table, 47, 79, 79, 1);
    const Color_RGBA darkslategrey(color_table, 47, 79, 79, 1);
    const Color_RGBA darkturquoise(color_table, 0, 206, 209, 1);
    const Color_RGBA darkviolet(color_table, 148, 0, 211, 1);
    const Color_RGBA deeppink(color_table, 255, 20, 147, 1);
    const Color_RGBA deepskyblue(color_table, 0, 191, 255, 1);
    const Color_RGBA dimgray(color_table, 105, 105, 105, 1);
    const Color_RGBA dimgrey(color_table, 105, 105, 105, 1);
    const Color_RGBA dodgerblue(color_table, 30, 144, 255, 1);
    const Color_RGBA firebrick(color_table, 178, 34, 34, 1);
    const Color_RGBA floralwhite(color_table, 255, 250, 240, 1);
    const Color_RGBA forestgreen(color_table, 34, 139, 34, 1);
    const Color_RGBA magenta(color_table, 255, 0, 255, 1);
    const Color_RGBA fuchsia(color_table, 255, 0, 255, 1);
    const Color_RGBA gainsboro(color_table, 220, 220, 220, 1);
    const Color_RGBA ghostwhite(color_table, 248, 248, 255, 1);
    const Color_RGBA gold(color_table, 255, 215, 0, 1);
    const Color_RGBA goldenrod(color_table, 218, 165, 32, 1);
    const Color_RGBA gray(color_table, 128, 128, 128, 1);
    const Color_RGBA grey(color_table, 128, 128, 128, 1);
    const Color_RGBA green(color_table, 0, 128, 0, 1);
    const Color_RGBA greenyellow(color_table, 173, 255, 47, 1);
    const Color_RGBA honeydew(color_table, 240, 255, 240, 1);
    const Color_RGBA hotpink(color_table, 255, 105, 180, 1);
    const Color_RGBA indianred(color_table, 205, 92, 92, 1);
    const Color_RGBA indigo(color_table, 75, 0, 130, 1);
    const Color_RGBA ivory(color_table, 255, 255, 240, 1);
    const Color_RGBA khaki(color_table, 240, 230, 140, 1);
    const Color_RGBA lavender(color_table, 230, 230, 250, 1);
    const Color_RGBA lavenderblush(color_table, 255, 240, 245, 1);
    const Color_RGBA lawngreen(color_table, 124, 252, 0, 1);
    const Color_RGBA lemonchiffon(color_table, 255, 250, 205, 1);
    const Color_RGBA lightblue(color_table, 173, 216, 230, 1);
    const Color_RGBA lightcoral(color_table, 240, 128, 128, 1);
    const Color_RGBA lightcyan(color_table, 224, 255, 255, 1);
    const Color_RGBA lightgoldenrodyellow(color_table, 250, 250, 210, 1);
    const Color_RGBA lightgray(color_table, 211, 211, 211, 1);
    const Color_RGBA lightgrey(color_table, 211, 211, 211, 1);
    const Color_RGBA lightgreen(color_table, 144, 238, 144, 1);
    const Color_RGBA lightpink(color_table, 255, 182, 193, 1);
    const Color_RGBA lightsalmon(color_table, 255, 160, 122, 1);
    const Color_RGBA lightseagreen(color_table, 32, 178, 170, 1);
    const Color_RGBA lightskyblue(color_table, 135, 206, 250, 1);
    const Color_RGBA lightslategray(color_table, 119, 136, 153, 1);
    const Color_RGBA lightslategrey(color_table, 119, 136, 153, 1);
    const Color_RGBA lightsteelblue(color_table, 176, 196, 222, 1);
    const Color_RGBA lightyellow(color_table, 255, 255, 224, 1);
    const Color_RGBA lime(color_table, 0, 255, 0, 1);
    const Color_RGBA limegreen(color_table, 50, 205, 50, 1);
    const Color_RGBA linen(color_table, 250, 240, 230, 1);
    const Color_RGBA maroon(color_table, 128, 0, 0, 1);
    const Color_RGBA mediumaquamarine(color_table, 102, 205, 170, 1);
    const Color_RGBA mediumblue(color_table, 0, 0, 205, 1);
    const Color_RGBA mediumorchid(color_table, 186, 85, 211, 1);
    const Color_RGBA mediumpurple(color_table, 147, 112, 219, 1);
    const Color_RGBA mediumseagreen(color_table, 60, 179, 113, 1);
    const Color_RGBA mediumslateblue(color_table, 123, 104, 238, 1);
    const Color_RGBA mediumspringgreen(color_table, 0, 250, 154, 1);
    const Color_RGBA mediumturquoise(color_table, 72, 209, 204, 1);
    const Color_RGBA mediumvioletred(color_table, 199, 21, 133, 1);
    const Color_RGBA midnightblue(color_table, 25, 25, 112, 1);
    const Color_RGBA mintcream(color_table, 245, 255, 250, 1);
    const Color_RGBA mistyrose(color_table, 255, 228, 225, 1);
    const Color_RGBA moccasin(color_table, 255, 228, 181, 1);
    const Color_RGBA navajowhite(color_table, 255, 222, 173, 1);
    const Color_RGBA navy(color_table, 0, 0, 128, 1);
    const Color_RGBA oldlace(color_table, 253, 245, 230, 1);
    const Color_RGBA olive(color_table, 128, 128, 0, 1);
    const Color_RGBA olivedrab(color_table, 107, 142, 35, 1);
    const Color_RGBA orange(color_table, 255, 165, 0, 1);
    const Color_RGBA orangered(color_table, 255, 69, 0, 1);
    const Color_RGBA orchid(color_table, 218, 112, 214, 1);
    const Color_RGBA palegoldenrod(color_table, 238, 232, 170, 1);
    const Color_RGBA palegreen(color_table, 152, 251, 152, 1);
    const Color_RGBA paleturquoise(color_table, 175, 238, 238, 1);
    const Color_RGBA palevioletred(color_table, 219, 112, 147, 1);
    const Color_RGBA papayawhip(color_table, 255, 239, 213, 1);
    const Color_RGBA peachpuff(color_table, 255, 218, 185, 1);
    const Color_RGBA peru(color_table, 205, 133, 63, 1);
    const Color_RGBA pink(color_table, 255, 192, 203, 1);
    const Color_RGBA plum(color_table, 221, 160, 221, 1);
    const Color_RGBA powderblue(color_table, 176, 224, 230, 1);
    const Color_RGBA purple(color_table, 128, 0, 128, 1);
    const Color_RGBA red(color_table, 255, 0, 0, 1);
    const Color_RGBA rosybrown(color_table, 188, 143, 143, 1);
    const Color_RGBA royalblue(color_table, 65, 105, 225, 1);
    const Color_RGBA saddlebrown(color_table, 139, 69, 19, 1);
    const Color_RGBA salmon(color_table, 250, 128, 114, 1);
    const Color_RGBA sandybrown(color_table, 244, 164, 96, 1);
    const Color_RGBA seagreen(color_table, 46, 139, 87, 1);
    const Color_RGBA seashell(color_table, 255, 245, 238, 1);
    const Color_RGBA sienna(color_table, 160, 82, 45, 1);
    const Color_RGBA silver(color_table, 192, 192, 192, 1);
    const Color_RGBA skyblue(color_table, 135, 206, 235, 1);
    const Color_RGBA slateblue(color_table, 106, 90, 205, 1);
    const Color_RGBA slategray(color_table, 112, 128, 144, 1);
    const Color_RGBA slategrey(color_table, 112, 128, 144, 1);
    const Color_RGBA snow(color_table, 255, 250, 250, 1);
    const Color_RGBA springgreen(color_table, 0, 255, 127, 1);
    const Color_RGBA steelblue(color_table, 70, 130, 180, 1);
    const Color_RGBA tan(color_table, 210, 180, 140, 1);
    const Color_RGBA teal(color_table, 0, 128, 128, 1);
    const Color_RGBA thistle(color_table, 216, 191, 216, 1);
    const Color_RGBA tomato(color_table, 255, 99, 71, 1);
    const Color_RGBA turquoise(color_table, 64, 224, 208, 1);
    const Color_RGBA violet(color_table, 238, 130, 238, 1);
    const Color_RGBA wheat(color_table, 245, 222, 179, 1);
    const Color_RGBA white(color_table, 255, 255, 255, 1);
    const Color_RGBA whitesmoke(color_table, 245, 245, 245, 1);
    const Color_RGBA yellow(color_table, 255, 255, 0, 1);
    const Color_RGBA yellowgreen(color_table, 154, 205, 50, 1);
    const Color_RGBA rebeccapurple(color_table, 102, 51, 153, 1);
    const Color_RGBA transparent(color_table, 0, 0, 0, 0);
  }

  const std::map<const int, const char*> colors_to_names {
    { 240 * 0x10000 + 248 * 0x100 + 255, ColorNames::aliceblue },
    { 250 * 0x10000 + 235 * 0x100 + 215, ColorNames::antiquewhite },
    {   0 * 0x10000 + 255 * 0x100 + 255, ColorNames::cyan },
    { 127 * 0x10000 + 255 * 0x100 + 212, ColorNames::aquamarine },
    { 240 * 0x10000 + 255 * 0x100 + 255, ColorNames::azure },
    { 245 * 0x10000 + 245 * 0x100 + 220, ColorNames::beige },
    { 255 * 0x10000 + 228 * 0x100 + 196, ColorNames::bisque },
    {   0 * 0x10000 +   0 * 0x100 +   0, ColorNames::black },
    { 255 * 0x10000 + 235 * 0x100 + 205, ColorNames::blanchedalmond },
    {   0 * 0x10000 +   0 * 0x100 + 255, ColorNames::blue },
    { 138 * 0x10000 +  43 * 0x100 + 226, ColorNames::blueviolet },
    { 165 * 0x10000 +  42 * 0x100 +  42, ColorNames::brown },
    { 222 * 0x10000 + 184 * 0x100 + 135, ColorNames::burlywood },
    {  95 * 0x10000 + 158 * 0x100 + 160, ColorNames::cadetblue },
    { 127 * 0x10000 + 255 * 0x100 +   0, ColorNames::chartreuse },
    { 210 * 0x10000 + 105 * 0x100 +  30, ColorNames::chocolate },
    { 255 * 0x10000 + 127 * 0x100 +  80, ColorNames::coral },
    { 100 * 0x10000 + 149 * 0x100 + 237, ColorNames::cornflowerblue },
    { 255 * 0x10000 + 248 * 0x100 + 220, ColorNames::cornsilk },
    { 220 * 0x10000 +  20 * 0x100 +  60, ColorNames::crimson },
    {   0 * 0x10000 +   0 * 0x100 + 139, ColorNames::darkblue },
    {   0 * 0x10000 + 139 * 0x100 + 139, ColorNames::darkcyan },
    { 184 * 0x10000 + 134 * 0x100 +  11, ColorNames::darkgoldenrod },
    { 169 * 0x10000 + 169 * 0x100 + 169, ColorNames::darkgray },
    {   0 * 0x10000 + 100 * 0x100 +   0, ColorNames::darkgreen },
    { 189 * 0x10000 + 183 * 0x100 + 107, ColorNames::darkkhaki },
    { 139 * 0x10000 +   0 * 0x100 + 139, ColorNames::darkmagenta },
    {  85 * 0x10000 + 107 * 0x100 +  47, ColorNames::darkolivegreen },
    { 255 * 0x10000 + 140 * 0x100 +   0, ColorNames::darkorange },
    { 153 * 0x10000 +  50 * 0x100 + 204, ColorNames::darkorchid },
    { 139 * 0x10000 +   0 * 0x100 +   0, ColorNames::darkred },
    { 233 * 0x10000 + 150 * 0x100 + 122, ColorNames::darksalmon },
    { 143 * 0x10000 + 188 * 0x100 + 143, ColorNames::darkseagreen },
    {  72 * 0x10000 +  61 * 0x100 + 139, ColorNames::darkslateblue },
    {  47 * 0x10000 +  79 * 0x100 +  79, ColorNames::darkslategray },
    {   0 * 0x10000 + 206 * 0x100 + 209, ColorNames::darkturquoise },
    { 148 * 0x10000 +   0 * 0x100 + 211, ColorNames::darkviolet },
    { 255 * 0x10000 +  20 * 0x100 + 147, ColorNames::deeppink },
    {   0 * 0x10000 + 191 * 0x100 + 255, ColorNames::deepskyblue },
    { 105 * 0x10000 + 105 * 0x100 + 105, ColorNames::dimgray },
    {  30 * 0x10000 + 144 * 0x100 + 255, ColorNames::dodgerblue },
    { 178 * 0x10000 +  34 * 0x100 +  34, ColorNames::firebrick },
    { 255 * 0x10000 + 250 * 0x100 + 240, ColorNames::floralwhite },
    {  34 * 0x10000 + 139 * 0x100 +  34, ColorNames::forestgreen },
    { 255 * 0x10000 +   0 * 0x100 + 255, ColorNames::magenta },
    { 220 * 0x10000 + 220 * 0x100 + 220, ColorNames::gainsboro },
    { 248 * 0x10000 + 248 * 0x100 + 255, ColorNames::ghostwhite },
    { 255 * 0x10000 + 215 * 0x100 +   0, ColorNames::gold },
    { 218 * 0x10000 + 165 * 0x100 +  32, ColorNames::goldenrod },
    { 128 * 0x10000 + 128 * 0x100 + 128, ColorNames::gray },
    {   0 * 0x10000 + 128 * 0x100 +   0, ColorNames::green },
    { 173 * 0x10000 + 255 * 0x100 +  47, ColorNames::greenyellow },
    { 240 * 0x10000 + 255 * 0x100 + 240, ColorNames::honeydew },
    { 255 * 0x10000 + 105 * 0x100 + 180, ColorNames::hotpink },
    { 205 * 0x10000 +  92 * 0x100 +  92, ColorNames::indianred },
    {  75 * 0x10000 +   0 * 0x100 + 130, ColorNames::indigo },
    { 255 * 0x10000 + 255 * 0x100 + 240, ColorNames::ivory },
    { 240 * 0x10000 + 230 * 0x100 + 140, ColorNames::khaki },
    { 230 * 0x10000 + 230 * 0x100 + 250, ColorNames::lavender },
    { 255 * 0x10000 + 240 * 0x100 + 245, ColorNames::lavenderblush },
    { 124 * 0x10000 + 252 * 0x100 +   0, ColorNames::lawngreen },
    { 255 * 0x10000 + 250 * 0x100 + 205, ColorNames::lemonchiffon },
    { 173 * 0x10000 + 216 * 0x100 + 230, ColorNames::lightblue },
    { 240 * 0x10000 + 128 * 0x100 + 128, ColorNames::lightcoral },
    { 224 * 0x10000 + 255 * 0x100 + 255, ColorNames::lightcyan },
    { 250 * 0x10000 + 250 * 0x100 + 210, ColorNames::lightgoldenrodyellow },
    { 211 * 0x10000 + 211 * 0x100 + 211, ColorNames::lightgray },
    { 144 * 0x10000 + 238 * 0x100 + 144, ColorNames::lightgreen },
    { 255 * 0x10000 + 182 * 0x100 + 193, ColorNames::lightpink },
    { 255 * 0x10000 + 160 * 0x100 + 122, ColorNames::lightsalmon },
    {  32 * 0x10000 + 178 * 0x100 + 170, ColorNames::lightseagreen },
    { 135 * 0x10000 + 206 * 0x100 + 250, ColorNames::lightskyblue },
    { 119 * 0x10000 + 136 * 0x100 + 153, ColorNames::lightslategray },
    { 176 * 0x10000 + 196 * 0x100 + 222, ColorNames::lightsteelblue },
    { 255 * 0x10000 + 255 * 0x100 + 224, ColorNames::lightyellow },
    {   0 * 0x10000 + 255 * 0x100 +   0, ColorNames::lime },
    {  50 * 0x10000 + 205 * 0x100 +  50, ColorNames::limegreen },
    { 250 * 0x10000 + 240 * 0x100 + 230, ColorNames::linen },
    { 128 * 0x10000 +   0 * 0x100 +   0, ColorNames::maroon },
    { 102 * 0x10000 + 205 * 0x100 + 170, ColorNames::mediumaquamarine },
    {   0 * 0x10000 +   0 * 0x100 + 205, ColorNames::mediumblue },
    { 186 * 0x10000 +  85 * 0x100 + 211, ColorNames::mediumorchid },
    { 147 * 0x10000 + 112 * 0x100 + 219, ColorNames::mediumpurple },
    {  60 * 0x10000 + 179 * 0x100 + 113, ColorNames::mediumseagreen },
    { 123 * 0x10000 + 104 * 0x100 + 238, ColorNames::mediumslateblue },
    {   0 * 0x10000 + 250 * 0x100 + 154, ColorNames::mediumspringgreen },
    {  72 * 0x10000 + 209 * 0x100 + 204, ColorNames::mediumturquoise },
    { 199 * 0x10000 +  21 * 0x100 + 133, ColorNames::mediumvioletred },
    {  25 * 0x10000 +  25 * 0x100 + 112, ColorNames::midnightblue },
    { 245 * 0x10000 + 255 * 0x100 + 250, ColorNames::mintcream },
    { 255 * 0x10000 + 228 * 0x100 + 225, ColorNames::mistyrose },
    { 255 * 0x10000 + 228 * 0x100 + 181, ColorNames::moccasin },
    { 255 * 0x10000 + 222 * 0x100 + 173, ColorNames::navajowhite },
    {   0 * 0x10000 +   0 * 0x100 + 128, ColorNames::navy },
    { 253 * 0x10000 + 245 * 0x100 + 230, ColorNames::oldlace },
    { 128 * 0x10000 + 128 * 0x100 +   0, ColorNames::olive },
    { 107 * 0x10000 + 142 * 0x100 +  35, ColorNames::olivedrab },
    { 255 * 0x10000 + 165 * 0x100 +   0, ColorNames::orange },
    { 255 * 0x10000 +  69 * 0x100 +   0, ColorNames::orangered },
    { 218 * 0x10000 + 112 * 0x100 + 214, ColorNames::orchid },
    { 238 * 0x10000 + 232 * 0x100 + 170, ColorNames::palegoldenrod },
    { 152 * 0x10000 + 251 * 0x100 + 152, ColorNames::palegreen },
    { 175 * 0x10000 + 238 * 0x100 + 238, ColorNames::paleturquoise },
    { 219 * 0x10000 + 112 * 0x100 + 147, ColorNames::palevioletred },
    { 255 * 0x10000 + 239 * 0x100 + 213, ColorNames::papayawhip },
    { 255 * 0x10000 + 218 * 0x100 + 185, ColorNames::peachpuff },
    { 205 * 0x10000 + 133 * 0x100 +  63, ColorNames::peru },
    { 255 * 0x10000 + 192 * 0x100 + 203, ColorNames::pink },
    { 221 * 0x10000 + 160 * 0x100 + 221, ColorNames::plum },
    { 176 * 0x10000 + 224 * 0x100 + 230, ColorNames::powderblue },
    { 128 * 0x10000 +   0 * 0x100 + 128, ColorNames::purple },
    { 255 * 0x10000 +   0 * 0x100 +   0, ColorNames::red },
    { 188 * 0x10000 + 143 * 0x100 + 143, ColorNames::rosybrown },
    {  65 * 0x10000 + 105 * 0x100 + 225, ColorNames::royalblue },
    { 139 * 0x10000 +  69 * 0x100 +  19, ColorNames::saddlebrown },
    { 250 * 0x10000 + 128 * 0x100 + 114, ColorNames::salmon },
    { 244 * 0x10000 + 164 * 0x100 +  96, ColorNames::sandybrown },
    {  46 * 0x10000 + 139 * 0x100 +  87, ColorNames::seagreen },
    { 255 * 0x10000 + 245 * 0x100 + 238, ColorNames::seashell },
    { 160 * 0x10000 +  82 * 0x100 +  45, ColorNames::sienna },
    { 192 * 0x10000 + 192 * 0x100 + 192, ColorNames::silver },
    { 135 * 0x10000 + 206 * 0x100 + 235, ColorNames::skyblue },
    { 106 * 0x10000 +  90 * 0x100 + 205, ColorNames::slateblue },
    { 112 * 0x10000 + 128 * 0x100 + 144, ColorNames::slategray },
    { 255 * 0x10000 + 250 * 0x100 + 250, ColorNames::snow },
    {   0 * 0x10000 + 255 * 0x100 + 127, ColorNames::springgreen },
    {  70 * 0x10000 + 130 * 0x100 + 180, ColorNames::steelblue },
    { 210 * 0x10000 + 180 * 0x100 + 140, ColorNames::tan },
    {   0 * 0x10000 + 128 * 0x100 + 128, ColorNames::teal },
    { 216 * 0x10000 + 191 * 0x100 + 216, ColorNames::thistle },
    { 255 * 0x10000 +  99 * 0x100 +  71, ColorNames::tomato },
    {  64 * 0x10000 + 224 * 0x100 + 208, ColorNames::turquoise },
    { 238 * 0x10000 + 130 * 0x100 + 238, ColorNames::violet },
    { 245 * 0x10000 + 222 * 0x100 + 179, ColorNames::wheat },
    { 255 * 0x10000 + 255 * 0x100 + 255, ColorNames::white },
    { 245 * 0x10000 + 245 * 0x100 + 245, ColorNames::whitesmoke },
    { 255 * 0x10000 + 255 * 0x100 +   0, ColorNames::yellow },
    { 154 * 0x10000 + 205 * 0x100 +  50, ColorNames::yellowgreen },
    { 102 * 0x10000 +  51 * 0x100 + 153, ColorNames::rebeccapurple }
  };

  const std::map<const char*, Color_RGBA_Ptr_Const, map_cmp_str> names_to_colors
  {
    { ColorNames::aliceblue, &Colors::aliceblue },
    { ColorNames::antiquewhite, &Colors::antiquewhite },
    { ColorNames::cyan, &Colors::cyan },
    { ColorNames::aqua, &Colors::aqua },
    { ColorNames::aquamarine, &Colors::aquamarine },
    { ColorNames::azure, &Colors::azure },
    { ColorNames::beige, &Colors::beige },
    { ColorNames::bisque, &Colors::bisque },
    { ColorNames::black, &Colors::black },
    { ColorNames::blanchedalmond, &Colors::blanchedalmond },
    { ColorNames::blue, &Colors::blue },
    { ColorNames::blueviolet, &Colors::blueviolet },
    { ColorNames::brown, &Colors::brown },
    { ColorNames::burlywood, &Colors::burlywood },
    { ColorNames::cadetblue, &Colors::cadetblue },
    { ColorNames::chartreuse, &Colors::chartreuse },
    { ColorNames::chocolate, &Colors::chocolate },
    { ColorNames::coral, &Colors::coral },
    { ColorNames::cornflowerblue, &Colors::cornflowerblue },
    { ColorNames::cornsilk, &Colors::cornsilk },
    { ColorNames::crimson, &Colors::crimson },
    { ColorNames::darkblue, &Colors::darkblue },
    { ColorNames::darkcyan, &Colors::darkcyan },
    { ColorNames::darkgoldenrod, &Colors::darkgoldenrod },
    { ColorNames::darkgray, &Colors::darkgray },
    { ColorNames::darkgrey, &Colors::darkgrey },
    { ColorNames::darkgreen, &Colors::darkgreen },
    { ColorNames::darkkhaki, &Colors::darkkhaki },
    { ColorNames::darkmagenta, &Colors::darkmagenta },
    { ColorNames::darkolivegreen, &Colors::darkolivegreen },
    { ColorNames::darkorange, &Colors::darkorange },
    { ColorNames::darkorchid, &Colors::darkorchid },
    { ColorNames::darkred, &Colors::darkred },
    { ColorNames::darksalmon, &Colors::darksalmon },
    { ColorNames::darkseagreen, &Colors::darkseagreen },
    { ColorNames::darkslateblue, &Colors::darkslateblue },
    { ColorNames::darkslategray, &Colors::darkslategray },
    { ColorNames::darkslategrey, &Colors::darkslategrey },
    { ColorNames::darkturquoise, &Colors::darkturquoise },
    { ColorNames::darkviolet, &Colors::darkviolet },
    { ColorNames::deeppink, &Colors::deeppink },
    { ColorNames::deepskyblue, &Colors::deepskyblue },
    { ColorNames::dimgray, &Colors::dimgray },
    { ColorNames::dimgrey, &Colors::dimgrey },
    { ColorNames::dodgerblue, &Colors::dodgerblue },
    { ColorNames::firebrick, &Colors::firebrick },
    { ColorNames::floralwhite, &Colors::floralwhite },
    { ColorNames::forestgreen, &Colors::forestgreen },
    { ColorNames::magenta, &Colors::magenta },
    { ColorNames::fuchsia, &Colors::fuchsia },
    { ColorNames::gainsboro, &Colors::gainsboro },
    { ColorNames::ghostwhite, &Colors::ghostwhite },
    { ColorNames::gold, &Colors::gold },
    { ColorNames::goldenrod, &Colors::goldenrod },
    { ColorNames::gray, &Colors::gray },
    { ColorNames::grey, &Colors::grey },
    { ColorNames::green, &Colors::green },
    { ColorNames::greenyellow, &Colors::greenyellow },
    { ColorNames::honeydew, &Colors::honeydew },
    { ColorNames::hotpink, &Colors::hotpink },
    { ColorNames::indianred, &Colors::indianred },
    { ColorNames::indigo, &Colors::indigo },
    { ColorNames::ivory, &Colors::ivory },
    { ColorNames::khaki, &Colors::khaki },
    { ColorNames::lavender, &Colors::lavender },
    { ColorNames::lavenderblush, &Colors::lavenderblush },
    { ColorNames::lawngreen, &Colors::lawngreen },
    { ColorNames::lemonchiffon, &Colors::lemonchiffon },
    { ColorNames::lightblue, &Colors::lightblue },
    { ColorNames::lightcoral, &Colors::lightcoral },
    { ColorNames::lightcyan, &Colors::lightcyan },
    { ColorNames::lightgoldenrodyellow, &Colors::lightgoldenrodyellow },
    { ColorNames::lightgray, &Colors::lightgray },
    { ColorNames::lightgrey, &Colors::lightgrey },
    { ColorNames::lightgreen, &Colors::lightgreen },
    { ColorNames::lightpink, &Colors::lightpink },
    { ColorNames::lightsalmon, &Colors::lightsalmon },
    { ColorNames::lightseagreen, &Colors::lightseagreen },
    { ColorNames::lightskyblue, &Colors::lightskyblue },
    { ColorNames::lightslategray, &Colors::lightslategray },
    { ColorNames::lightslategrey, &Colors::lightslategrey },
    { ColorNames::lightsteelblue, &Colors::lightsteelblue },
    { ColorNames::lightyellow, &Colors::lightyellow },
    { ColorNames::lime, &Colors::lime },
    { ColorNames::limegreen, &Colors::limegreen },
    { ColorNames::linen, &Colors::linen },
    { ColorNames::maroon, &Colors::maroon },
    { ColorNames::mediumaquamarine, &Colors::mediumaquamarine },
    { ColorNames::mediumblue, &Colors::mediumblue },
    { ColorNames::mediumorchid, &Colors::mediumorchid },
    { ColorNames::mediumpurple, &Colors::mediumpurple },
    { ColorNames::mediumseagreen, &Colors::mediumseagreen },
    { ColorNames::mediumslateblue, &Colors::mediumslateblue },
    { ColorNames::mediumspringgreen, &Colors::mediumspringgreen },
    { ColorNames::mediumturquoise, &Colors::mediumturquoise },
    { ColorNames::mediumvioletred, &Colors::mediumvioletred },
    { ColorNames::midnightblue, &Colors::midnightblue },
    { ColorNames::mintcream, &Colors::mintcream },
    { ColorNames::mistyrose, &Colors::mistyrose },
    { ColorNames::moccasin, &Colors::moccasin },
    { ColorNames::navajowhite, &Colors::navajowhite },
    { ColorNames::navy, &Colors::navy },
    { ColorNames::oldlace, &Colors::oldlace },
    { ColorNames::olive, &Colors::olive },
    { ColorNames::olivedrab, &Colors::olivedrab },
    { ColorNames::orange, &Colors::orange },
    { ColorNames::orangered, &Colors::orangered },
    { ColorNames::orchid, &Colors::orchid },
    { ColorNames::palegoldenrod, &Colors::palegoldenrod },
    { ColorNames::palegreen, &Colors::palegreen },
    { ColorNames::paleturquoise, &Colors::paleturquoise },
    { ColorNames::palevioletred, &Colors::palevioletred },
    { ColorNames::papayawhip, &Colors::papayawhip },
    { ColorNames::peachpuff, &Colors::peachpuff },
    { ColorNames::peru, &Colors::peru },
    { ColorNames::pink, &Colors::pink },
    { ColorNames::plum, &Colors::plum },
    { ColorNames::powderblue, &Colors::powderblue },
    { ColorNames::purple, &Colors::purple },
    { ColorNames::red, &Colors::red },
    { ColorNames::rosybrown, &Colors::rosybrown },
    { ColorNames::royalblue, &Colors::royalblue },
    { ColorNames::saddlebrown, &Colors::saddlebrown },
    { ColorNames::salmon, &Colors::salmon },
    { ColorNames::sandybrown, &Colors::sandybrown },
    { ColorNames::seagreen, &Colors::seagreen },
    { ColorNames::seashell, &Colors::seashell },
    { ColorNames::sienna, &Colors::sienna },
    { ColorNames::silver, &Colors::silver },
    { ColorNames::skyblue, &Colors::skyblue },
    { ColorNames::slateblue, &Colors::slateblue },
    { ColorNames::slategray, &Colors::slategray },
    { ColorNames::slategrey, &Colors::slategrey },
    { ColorNames::snow, &Colors::snow },
    { ColorNames::springgreen, &Colors::springgreen },
    { ColorNames::steelblue, &Colors::steelblue },
    { ColorNames::tan, &Colors::tan },
    { ColorNames::teal, &Colors::teal },
    { ColorNames::thistle, &Colors::thistle },
    { ColorNames::tomato, &Colors::tomato },
    { ColorNames::turquoise, &Colors::turquoise },
    { ColorNames::violet, &Colors::violet },
    { ColorNames::wheat, &Colors::wheat },
    { ColorNames::white, &Colors::white },
    { ColorNames::whitesmoke, &Colors::whitesmoke },
    { ColorNames::yellow, &Colors::yellow },
    { ColorNames::yellowgreen, &Colors::yellowgreen },
    { ColorNames::rebeccapurple, &Colors::rebeccapurple },
    { ColorNames::transparent, &Colors::transparent }
  };

  Color_RGBA_Ptr_Const name_to_color(const char* key)
  {
    return name_to_color(std::string(key));
  }

  Color_RGBA_Ptr_Const name_to_color(const std::string& key)
  {
    // case insensitive lookup.  See #2462
    std::string lower{key};
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    auto p = names_to_colors.find(lower.c_str());
    if (p != names_to_colors.end()) {
      return p->second;
    }
    return 0;
  }

  const char* color_to_name(const int key)
  {
    auto p = colors_to_names.find(key);
    if (p != colors_to_names.end()) {
      return p->second;
    }
    return 0;
  }

  const char* color_to_name(const double key)
  {
    return color_to_name((int)key);
  }

  const char* color_to_name(const Color_RGBA& c)
  {
    double key = c.r() * 0x10000
               + c.g() * 0x100
               + c.b();
    return color_to_name(key);
  }

}
/*
cencoder.c - c source to a base64 encoding algorithm implementation

This is part of the libb64 project, and has been placed in the public domain.
For details, see http://sourceforge.net/projects/libb64
*/

namespace base64 {

void base64_init_encodestate(base64_encodestate* state_in)
{
	state_in->step = step_A;
	state_in->result = 0;
	state_in->stepcount = 0;
}

char base64_encode_value(char value_in)
{
	static const char* encoding = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	if (value_in > 63) return '=';
	return encoding[(int)value_in];
}

int base64_encode_block(const char* plaintext_in, int length_in, char* code_out, base64_encodestate* state_in)
{
	const char* plainchar = plaintext_in;
	const char* const plaintextend = plaintext_in + length_in;
	char* codechar = code_out;
	char result;
	char fragment;

	result = state_in->result;

	switch (state_in->step)
	{
		while (1)
		{
	case step_A:
			if (plainchar == plaintextend)
			{
				state_in->result = result;
				state_in->step = step_A;
				return (int)(codechar - code_out);
			}
			fragment = *plainchar++;
			result = (fragment & 0x0fc) >> 2;
			*codechar++ = base64_encode_value(result);
			result = (fragment & 0x003) << 4;
			#ifndef _MSC_VER
				/* fall through */
			#endif
	case step_B:
			if (plainchar == plaintextend)
			{
				state_in->result = result;
				state_in->step = step_B;
				return (int)(codechar - code_out);
			}
			fragment = *plainchar++;
			result |= (fragment & 0x0f0) >> 4;
			*codechar++ = base64_encode_value(result);
			result = (fragment & 0x00f) << 2;
			#ifndef _MSC_VER
				/* fall through */
			#endif
	case step_C:
			if (plainchar == plaintextend)
			{
				state_in->result = result;
				state_in->step = step_C;
				return (int)(codechar - code_out);
			}
			fragment = *plainchar++;
			result |= (fragment & 0x0c0) >> 6;
			*codechar++ = base64_encode_value(result);
			result  = (fragment & 0x03f) >> 0;
			*codechar++ = base64_encode_value(result);

			++(state_in->stepcount);
		}
	}
	/* control should not reach here */
	return (int)(codechar - code_out);
}

int base64_encode_blockend(char* code_out, base64_encodestate* state_in)
{
	char* codechar = code_out;

	switch (state_in->step)
	{
	case step_B:
		*codechar++ = base64_encode_value(state_in->result);
		*codechar++ = '=';
		*codechar++ = '=';
		break;
	case step_C:
		*codechar++ = base64_encode_value(state_in->result);
		*codechar++ = '=';
		break;
	case step_A:
		break;
	}
	*codechar++ = '\n';

	return (int)(codechar - code_out);
}

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  Output::Output(Sass_Output_Options& opt)
  : Inspect(Emitter(opt)),
    charset(""),
    top_nodes(0)
  {}

  Output::~Output() { }

  void Output::fallback_impl(AST_Node_Ptr n)
  {
    return n->perform(this);
  }

  void Output::operator()(Number_Ptr n)
  {
    // check for a valid unit here
    // includes result for reporting
    if (!n->is_valid_css_unit()) {
      // should be handle in check_expression
      throw Exception::InvalidValue({}, *n);
    }
    // use values to_string facility
    std::string res = n->to_string(opt);
    // output the final token
    append_token(res, n);
  }

  void Output::operator()(Import_Ptr imp)
  {
    top_nodes.push_back(imp);
  }

  void Output::operator()(Map_Ptr m)
  {
    // should be handle in check_expression
    throw Exception::InvalidValue({}, *m);
  }

  OutputBuffer Output::get_buffer(void)
  {

    Emitter emitter(opt);
    Inspect inspect(emitter);

    size_t size_nodes = top_nodes.size();
    for (size_t i = 0; i < size_nodes; i++) {
      top_nodes[i]->perform(&inspect);
      inspect.append_mandatory_linefeed();
    }

    // flush scheduled outputs
    // maybe omit semicolon if possible
    inspect.finalize(wbuf.buffer.size() == 0);
    // prepend buffer on top
    prepend_output(inspect.output());
    // make sure we end with a linefeed
    if (!ends_with(wbuf.buffer, opt.linefeed)) {
      // if the output is not completely empty
      if (!wbuf.buffer.empty()) append_string(opt.linefeed);
    }

    // search for unicode char
    for(const char& chr : wbuf.buffer) {
      // skip all ascii chars
      // static cast to unsigned to handle `char` being signed / unsigned
      if (static_cast<unsigned>(chr) < 128) continue;
      // declare the charset
      if (output_style() != COMPRESSED)
        charset = "@charset \"UTF-8\";"
                + std::string(opt.linefeed);
      else charset = "\xEF\xBB\xBF";
      // abort search
      break;
    }

    // add charset as first line, before comments and imports
    if (!charset.empty()) prepend_string(charset);

    return wbuf;

  }

  void Output::operator()(Comment_Ptr c)
  {
    std::string txt = c->text()->to_string(opt);
    // if (indentation && txt == "/**/") return;
    bool important = c->is_important();
    if (output_style() != COMPRESSED || important) {
      if (buffer().size() == 0) {
        top_nodes.push_back(c);
      } else {
        in_comment = true;
        append_indentation();
        c->text()->perform(this);
        in_comment = false;
        if (indentation == 0) {
          append_mandatory_linefeed();
        } else {
          append_optional_linefeed();
        }
      }
    }
  }

  void Output::operator()(Ruleset_Ptr r)
  {
    Selector_Obj s     = r->selector();
    Block_Obj    b     = r->block();

    // Filter out rulesets that aren't printable (process its children though)
    if (!Util::isPrintable(r, output_style())) {
      for (size_t i = 0, L = b->length(); i < L; ++i) {
        const Statement_Obj& stm = b->at(i);
        if (Cast<Has_Block>(stm)) {
          if (!Cast<Declaration>(stm)) {
            stm->perform(this);
          }
        }
      }
      return;
    }

    if (output_style() == NESTED) indentation += r->tabs();
    if (opt.source_comments) {
      std::stringstream ss;
      append_indentation();
      std::string path(File::abs2rel(r->pstate().path));
      ss << "/* line " << r->pstate().line + 1 << ", " << path << " */";
      append_string(ss.str());
      append_optional_linefeed();
    }
    scheduled_crutch = s;
    if (s) s->perform(this);
    append_scope_opener(b);
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj stm = b->at(i);
      bool bPrintExpression = true;
      // Check print conditions
      if (Declaration_Ptr dec = Cast<Declaration>(stm)) {
        if (String_Constant_Ptr valConst = Cast<String_Constant>(dec->value())) {
          std::string val(valConst->value());
          if (String_Quoted_Ptr qstr = Cast<String_Quoted>(valConst)) {
            if (!qstr->quote_mark() && val.empty()) {
              bPrintExpression = false;
            }
          }
        }
        else if (List_Ptr list = Cast<List>(dec->value())) {
          bool all_invisible = true;
          for (size_t list_i = 0, list_L = list->length(); list_i < list_L; ++list_i) {
            Expression_Ptr item = list->at(list_i);
            if (!item->is_invisible()) all_invisible = false;
          }
          if (all_invisible && !list->is_bracketed()) bPrintExpression = false;
        }
      }
      // Print if OK
      if (bPrintExpression) {
        stm->perform(this);
      }
    }
    if (output_style() == NESTED) indentation -= r->tabs();
    append_scope_closer(b);

  }
  void Output::operator()(Keyframe_Rule_Ptr r)
  {
    Block_Obj b = r->block();
    Selector_Obj v = r->name();

    if (!v.isNull()) {
      v->perform(this);
    }

    if (!b) {
      append_colon_separator();
      return;
    }

    append_scope_opener();
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj stm = b->at(i);
      stm->perform(this);
      if (i < L - 1) append_special_linefeed();
    }
    append_scope_closer();
  }

  void Output::operator()(Supports_Block_Ptr f)
  {
    if (f->is_invisible()) return;

    Supports_Condition_Obj c = f->condition();
    Block_Obj b              = f->block();

    // Filter out feature blocks that aren't printable (process its children though)
    if (!Util::isPrintable(f, output_style())) {
      for (size_t i = 0, L = b->length(); i < L; ++i) {
        Statement_Obj stm = b->at(i);
        if (Cast<Has_Block>(stm)) {
          stm->perform(this);
        }
      }
      return;
    }

    if (output_style() == NESTED) indentation += f->tabs();
    append_indentation();
    append_token("@supports", f);
    append_mandatory_space();
    c->perform(this);
    append_scope_opener();

    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj stm = b->at(i);
      stm->perform(this);
      if (i < L - 1) append_special_linefeed();
    }

    if (output_style() == NESTED) indentation -= f->tabs();

    append_scope_closer();

  }

  void Output::operator()(Media_Block_Ptr m)
  {
    if (m->is_invisible()) return;

    Block_Obj b     = m->block();

    // Filter out media blocks that aren't printable (process its children though)
    if (!Util::isPrintable(m, output_style())) {
      for (size_t i = 0, L = b->length(); i < L; ++i) {
        Statement_Obj stm = b->at(i);
        if (Cast<Has_Block>(stm)) {
          stm->perform(this);
        }
      }
      return;
    }
    if (output_style() == NESTED) indentation += m->tabs();
    append_indentation();
    append_token("@media", m);
    append_mandatory_space();
    in_media_block = true;
    m->media_queries()->perform(this);
    in_media_block = false;
    append_scope_opener();

    for (size_t i = 0, L = b->length(); i < L; ++i) {
      if (b->at(i)) {
      Statement_Obj stm = b->at(i);
        stm->perform(this);
      }
      if (i < L - 1) append_special_linefeed();
    }

    if (output_style() == NESTED) indentation -= m->tabs();
    append_scope_closer();
  }

  void Output::operator()(Directive_Ptr a)
  {
    std::string      kwd   = a->keyword();
    Selector_Obj   s     = a->selector();
    Expression_Obj v     = a->value();
    Block_Obj      b     = a->block();

    append_indentation();
    append_token(kwd, a);
    if (s) {
      append_mandatory_space();
      in_wrapped = true;
      s->perform(this);
      in_wrapped = false;
    }
    if (v) {
      append_mandatory_space();
      // ruby sass bug? should use options?
      append_token(v->to_string(/* opt */), v);
    }
    if (!b) {
      append_delimiter();
      return;
    }

    if (b->is_invisible() || b->length() == 0) {
      append_optional_space();
      return append_string("{}");
    }

    append_scope_opener();

    bool format = kwd != "@font-face";;

    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj stm = b->at(i);
      stm->perform(this);
      if (i < L - 1 && format) append_special_linefeed();
    }

    append_scope_closer();
  }

  void Output::operator()(String_Quoted_Ptr s)
  {
    if (s->quote_mark()) {
      append_token(quote(s->value(), s->quote_mark()), s);
    } else if (!in_comment) {
      append_token(string_to_output(s->value()), s);
    } else {
      append_token(s->value(), s);
    }
  }

  void Output::operator()(String_Constant_Ptr s)
  {
    std::string value(s->value());
    if (s->can_compress_whitespace() && output_style() == COMPRESSED) {
      value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());
    }
    if (!in_comment && !in_custom_property) {
      append_token(string_to_output(value), s);
    } else {
      append_token(value, s);
    }
  }

}

namespace Sass {

  template <typename T>
  Environment<T>::Environment(bool is_shadow)
  : local_frame_(environment_map<std::string, T>()),
    parent_(0), is_shadow_(false)
  { }
  template <typename T>
  Environment<T>::Environment(Environment<T>* env, bool is_shadow)
  : local_frame_(environment_map<std::string, T>()),
    parent_(env), is_shadow_(is_shadow)
  { }
  template <typename T>
  Environment<T>::Environment(Environment<T>& env, bool is_shadow)
  : local_frame_(environment_map<std::string, T>()),
    parent_(&env), is_shadow_(is_shadow)
  { }

  // link parent to create a stack
  template <typename T>
  void Environment<T>::link(Environment& env) { parent_ = &env; }
  template <typename T>
  void Environment<T>::link(Environment* env) { parent_ = env; }

  // this is used to find the global frame
  // which is the second last on the stack
  template <typename T>
  bool Environment<T>::is_lexical() const
  {
    return !! parent_ && parent_->parent_;
  }

  // only match the real root scope
  // there is still a parent around
  // not sure what it is actually use for
  // I guess we store functions etc. there
  template <typename T>
  bool Environment<T>::is_global() const
  {
    return parent_ && ! parent_->parent_;
  }

  template <typename T>
  environment_map<std::string, T>& Environment<T>::local_frame() {
    return local_frame_;
  }

  template <typename T>
  bool Environment<T>::has_local(const std::string& key) const
  { return local_frame_.find(key) != local_frame_.end(); }

  template <typename T> EnvResult
  Environment<T>::find_local(const std::string& key)
  {
    auto end = local_frame_.end();
    auto it = local_frame_.find(key);
    return EnvResult(it, it != end);
  }

  template <typename T>
  T& Environment<T>::get_local(const std::string& key)
  { return local_frame_[key]; }

  template <typename T>
  void Environment<T>::set_local(const std::string& key, const T& val)
  {
    local_frame_[key] = val;
  }
  template <typename T>
  void Environment<T>::set_local(const std::string& key, T&& val)
  {
    local_frame_[key] = val;
  }

  template <typename T>
  void Environment<T>::del_local(const std::string& key)
  { local_frame_.erase(key); }

  template <typename T>
  Environment<T>* Environment<T>::global_env()
  {
    Environment* cur = this;
    while (cur->is_lexical()) {
      cur = cur->parent_;
    }
    return cur;
  }

  template <typename T>
  bool Environment<T>::has_global(const std::string& key)
  { return global_env()->has(key); }

  template <typename T>
  T& Environment<T>::get_global(const std::string& key)
  { return (*global_env())[key]; }

  template <typename T>
  void Environment<T>::set_global(const std::string& key, const T& val)
  {
    global_env()->local_frame_[key] = val;
  }
  template <typename T>
  void Environment<T>::set_global(const std::string& key, T&& val)
  {
    global_env()->local_frame_[key] = val;
  }

  template <typename T>
  void Environment<T>::del_global(const std::string& key)
  { global_env()->local_frame_.erase(key); }

  template <typename T>
  Environment<T>* Environment<T>::lexical_env(const std::string& key)
  {
    Environment* cur = this;
    while (cur) {
      if (cur->has_local(key)) {
        return cur;
      }
      cur = cur->parent_;
    }
    return this;
  }

  // see if we have a lexical variable
  // move down the stack but stop before we
  // reach the global frame (is not included)
  template <typename T>
  bool Environment<T>::has_lexical(const std::string& key) const
  {
    auto cur = this;
    while (cur->is_lexical()) {
      if (cur->has_local(key)) return true;
      cur = cur->parent_;
    }
    return false;
  }

  // see if we have a lexical we could update
  // either update already existing lexical value
  // or if flag is set, we create one if no lexical found
  template <typename T>
  void Environment<T>::set_lexical(const std::string& key, const T& val)
  {
    Environment<T>* cur = this;
    bool shadow = false;
    while ((cur && cur->is_lexical()) || shadow) {
      EnvResult rv(cur->find_local(key));
      if (rv.found) {
        rv.it->second = val;
        return;
      }
      shadow = cur->is_shadow();
      cur = cur->parent_;
    }
    set_local(key, val);
  }
  // this one moves the value
  template <typename T>
  void Environment<T>::set_lexical(const std::string& key, T&& val)
  {
    Environment<T>* cur = this;
    bool shadow = false;
    while ((cur && cur->is_lexical()) || shadow) {
      EnvResult rv(cur->find_local(key));
      if (rv.found) {
        rv.it->second = val;
        return;
      }
      shadow = cur->is_shadow();
      cur = cur->parent_;
    }
    set_local(key, val);
  }

  // look on the full stack for key
  // include all scopes available
  template <typename T>
  bool Environment<T>::has(const std::string& key) const
  {
    auto cur = this;
    while (cur) {
      if (cur->has_local(key)) {
        return true;
      }
      cur = cur->parent_;
    }
    return false;
  }

  // look on the full stack for key
  // include all scopes available
  template <typename T> EnvResult
  Environment<T>::find(const std::string& key)
  {
    auto cur = this;
    while (true) {
      EnvResult rv(cur->find_local(key));
      if (rv.found) return rv;
      cur = cur->parent_;
      if (!cur) return rv;
    }
  };

  // use array access for getter and setter functions
  template <typename T>
  T& Environment<T>::get(const std::string& key)
  {
    auto cur = this;
    while (cur) {
      if (cur->has_local(key)) {
        return cur->get_local(key);
      }
      cur = cur->parent_;
    }
    return get_local(key);
  }

  // use array access for getter and setter functions
  template <typename T>
  T& Environment<T>::operator[](const std::string& key)
  {
    auto cur = this;
    while (cur) {
      if (cur->has_local(key)) {
        return cur->get_local(key);
      }
      cur = cur->parent_;
    }
    return get_local(key);
  }
/*
  #ifdef DEBUG
  template <typename T>
  size_t Environment<T>::print(std::string prefix)
  {
    size_t indent = 0;
    if (parent_) indent = parent_->print(prefix) + 1;
    std::cerr << prefix << std::string(indent, ' ') << "== " << this << std::endl;
    for (typename environment_map<std::string, T>::iterator i = local_frame_.begin(); i != local_frame_.end(); ++i) {
      if (!ends_with(i->first, "[f]") && !ends_with(i->first, "[f]4") && !ends_with(i->first, "[f]2")) {
        std::cerr << prefix << std::string(indent, ' ') << i->first << " " << i->second;
        if (Value_Ptr val = Cast<Value>(i->second))
        { std::cerr << " : " << val->to_string(); }
        std::cerr << std::endl;
      }
    }
    return indent ;
  }
  #endif
*/
  // compile implementation for AST_Node
  template class Environment<AST_Node_Obj>;

}


namespace Sass {

  namespace Functions {

    /////////////////
    // LIST FUNCTIONS
    /////////////////

    Signature keywords_sig = "keywords($args)";
    BUILT_IN(keywords)
    {
      List_Obj arglist = SASS_MEMORY_COPY(ARG("$args", List)); // copy
      Map_Obj result = SASS_MEMORY_NEW(Map, pstate, 1);
      for (size_t i = arglist->size(), L = arglist->length(); i < L; ++i) {
        Expression_Obj obj = arglist->at(i);
        Argument_Obj arg = (Argument_Ptr) obj.ptr(); // XXX
        std::string name = std::string(arg->name());
        name = name.erase(0, 1); // sanitize name (remove dollar sign)
        *result << std::make_pair(SASS_MEMORY_NEW(String_Quoted,
                 pstate, name),
                 arg->value());
      }
      return result.detach();
    }

    Signature length_sig = "length($list)";
    BUILT_IN(length)
    {
      if (Selector_List_Ptr sl = Cast<Selector_List>(env["$list"])) {
        return SASS_MEMORY_NEW(Number, pstate, (double)sl->length());
      }
      Expression_Ptr v = ARG("$list", Expression);
      if (v->concrete_type() == Expression::MAP) {
        Map_Ptr map = Cast<Map>(env["$list"]);
        return SASS_MEMORY_NEW(Number, pstate, (double)(map ? map->length() : 1));
      }
      if (v->concrete_type() == Expression::SELECTOR) {
        if (Compound_Selector_Ptr h = Cast<Compound_Selector>(v)) {
          return SASS_MEMORY_NEW(Number, pstate, (double)h->length());
        } else if (Selector_List_Ptr ls = Cast<Selector_List>(v)) {
          return SASS_MEMORY_NEW(Number, pstate, (double)ls->length());
        } else {
          return SASS_MEMORY_NEW(Number, pstate, 1);
        }
      }

      List_Ptr list = Cast<List>(env["$list"]);
      return SASS_MEMORY_NEW(Number,
                             pstate,
                             (double)(list ? list->size() : 1));
    }

    Signature nth_sig = "nth($list, $n)";
    BUILT_IN(nth)
    {
      double nr = ARGVAL("$n");
      Map_Ptr m = Cast<Map>(env["$list"]);
      if (Selector_List_Ptr sl = Cast<Selector_List>(env["$list"])) {
        size_t len = m ? m->length() : sl->length();
        bool empty = m ? m->empty() : sl->empty();
        if (empty) error("argument `$list` of `" + std::string(sig) + "` must not be empty", pstate, traces);
        double index = std::floor(nr < 0 ? len + nr : nr - 1);
        if (index < 0 || index > len - 1) error("index out of bounds for `" + std::string(sig) + "`", pstate, traces);
        // return (*sl)[static_cast<int>(index)];
        Listize listize;
        return Cast<Value>((*sl)[static_cast<int>(index)]->perform(&listize));
      }
      List_Obj l = Cast<List>(env["$list"]);
      if (nr == 0) error("argument `$n` of `" + std::string(sig) + "` must be non-zero", pstate, traces);
      // if the argument isn't a list, then wrap it in a singleton list
      if (!m && !l) {
        l = SASS_MEMORY_NEW(List, pstate, 1);
        l->append(ARG("$list", Expression));
      }
      size_t len = m ? m->length() : l->length();
      bool empty = m ? m->empty() : l->empty();
      if (empty) error("argument `$list` of `" + std::string(sig) + "` must not be empty", pstate, traces);
      double index = std::floor(nr < 0 ? len + nr : nr - 1);
      if (index < 0 || index > len - 1) error("index out of bounds for `" + std::string(sig) + "`", pstate, traces);

      if (m) {
        l = SASS_MEMORY_NEW(List, pstate, 2);
        l->append(m->keys()[static_cast<unsigned int>(index)]);
        l->append(m->at(m->keys()[static_cast<unsigned int>(index)]));
        return l.detach();
      }
      else {
        Value_Obj rv = l->value_at_index(static_cast<int>(index));
        rv->set_delayed(false);
        return rv.detach();
      }
    }

    Signature set_nth_sig = "set-nth($list, $n, $value)";
    BUILT_IN(set_nth)
    {
      Map_Obj m = Cast<Map>(env["$list"]);
      List_Obj l = Cast<List>(env["$list"]);
      Number_Obj n = ARG("$n", Number);
      Expression_Obj v = ARG("$value", Expression);
      if (!l) {
        l = SASS_MEMORY_NEW(List, pstate, 1);
        l->append(ARG("$list", Expression));
      }
      if (m) {
        l = m->to_list(pstate);
      }
      if (l->empty()) error("argument `$list` of `" + std::string(sig) + "` must not be empty", pstate, traces);
      double index = std::floor(n->value() < 0 ? l->length() + n->value() : n->value() - 1);
      if (index < 0 || index > l->length() - 1) error("index out of bounds for `" + std::string(sig) + "`", pstate, traces);
      List_Ptr result = SASS_MEMORY_NEW(List, pstate, l->length(), l->separator(), false, l->is_bracketed());
      for (size_t i = 0, L = l->length(); i < L; ++i) {
        result->append(((i == index) ? v : (*l)[i]));
      }
      return result;
    }

    Signature index_sig = "index($list, $value)";
    BUILT_IN(index)
    {
      Map_Obj m = Cast<Map>(env["$list"]);
      List_Obj l = Cast<List>(env["$list"]);
      Expression_Obj v = ARG("$value", Expression);
      if (!l) {
        l = SASS_MEMORY_NEW(List, pstate, 1);
        l->append(ARG("$list", Expression));
      }
      if (m) {
        l = m->to_list(pstate);
      }
      for (size_t i = 0, L = l->length(); i < L; ++i) {
        if (Operators::eq(l->value_at_index(i), v)) return SASS_MEMORY_NEW(Number, pstate, (double)(i+1));
      }
      return SASS_MEMORY_NEW(Null, pstate);
    }

    Signature join_sig = "join($list1, $list2, $separator: auto, $bracketed: auto)";
    BUILT_IN(join)
    {
      Map_Obj m1 = Cast<Map>(env["$list1"]);
      Map_Obj m2 = Cast<Map>(env["$list2"]);
      List_Obj l1 = Cast<List>(env["$list1"]);
      List_Obj l2 = Cast<List>(env["$list2"]);
      String_Constant_Obj sep = ARG("$separator", String_Constant);
      enum Sass_Separator sep_val = (l1 ? l1->separator() : SASS_SPACE);
      Value* bracketed = ARG("$bracketed", Value);
      bool is_bracketed = (l1 ? l1->is_bracketed() : false);
      if (!l1) {
        l1 = SASS_MEMORY_NEW(List, pstate, 1);
        l1->append(ARG("$list1", Expression));
        sep_val = (l2 ? l2->separator() : SASS_SPACE);
        is_bracketed = (l2 ? l2->is_bracketed() : false);
      }
      if (!l2) {
        l2 = SASS_MEMORY_NEW(List, pstate, 1);
        l2->append(ARG("$list2", Expression));
      }
      if (m1) {
        l1 = m1->to_list(pstate);
        sep_val = SASS_COMMA;
      }
      if (m2) {
        l2 = m2->to_list(pstate);
      }
      size_t len = l1->length() + l2->length();
      std::string sep_str = unquote(sep->value());
      if (sep_str == "space") sep_val = SASS_SPACE;
      else if (sep_str == "comma") sep_val = SASS_COMMA;
      else if (sep_str != "auto") error("argument `$separator` of `" + std::string(sig) + "` must be `space`, `comma`, or `auto`", pstate, traces);
      String_Constant_Obj bracketed_as_str = Cast<String_Constant>(bracketed);
      bool bracketed_is_auto = bracketed_as_str && unquote(bracketed_as_str->value()) == "auto";
      if (!bracketed_is_auto) {
        is_bracketed = !bracketed->is_false();
      }
      List_Obj result = SASS_MEMORY_NEW(List, pstate, len, sep_val, false, is_bracketed);
      result->concat(l1);
      result->concat(l2);
      return result.detach();
    }

    Signature append_sig = "append($list, $val, $separator: auto)";
    BUILT_IN(append)
    {
      Map_Obj m = Cast<Map>(env["$list"]);
      List_Obj l = Cast<List>(env["$list"]);
      Expression_Obj v = ARG("$val", Expression);
      if (Selector_List_Ptr sl = Cast<Selector_List>(env["$list"])) {
        Listize listize;
        l = Cast<List>(sl->perform(&listize));
      }
      String_Constant_Obj sep = ARG("$separator", String_Constant);
      if (!l) {
        l = SASS_MEMORY_NEW(List, pstate, 1);
        l->append(ARG("$list", Expression));
      }
      if (m) {
        l = m->to_list(pstate);
      }
      List_Ptr result = SASS_MEMORY_COPY(l);
      std::string sep_str(unquote(sep->value()));
      if (sep_str != "auto") { // check default first
        if (sep_str == "space") result->separator(SASS_SPACE);
        else if (sep_str == "comma") result->separator(SASS_COMMA);
        else error("argument `$separator` of `" + std::string(sig) + "` must be `space`, `comma`, or `auto`", pstate, traces);
      }
      if (l->is_arglist()) {
        result->append(SASS_MEMORY_NEW(Argument,
                                       v->pstate(),
                                       v,
                                       "",
                                       false,
                                       false));

      } else {
        result->append(v);
      }
      return result;
    }

    Signature zip_sig = "zip($lists...)";
    BUILT_IN(zip)
    {
      List_Obj arglist = SASS_MEMORY_COPY(ARG("$lists", List));
      size_t shortest = 0;
      for (size_t i = 0, L = arglist->length(); i < L; ++i) {
        List_Obj ith = Cast<List>(arglist->value_at_index(i));
        Map_Obj mith = Cast<Map>(arglist->value_at_index(i));
        if (!ith) {
          if (mith) {
            ith = mith->to_list(pstate);
          } else {
            ith = SASS_MEMORY_NEW(List, pstate, 1);
            ith->append(arglist->value_at_index(i));
          }
          if (arglist->is_arglist()) {
            Argument_Obj arg = (Argument_Ptr)(arglist->at(i).ptr()); // XXX
            arg->value(ith);
          } else {
            (*arglist)[i] = ith;
          }
        }
        shortest = (i ? std::min(shortest, ith->length()) : ith->length());
      }
      List_Ptr zippers = SASS_MEMORY_NEW(List, pstate, shortest, SASS_COMMA);
      size_t L = arglist->length();
      for (size_t i = 0; i < shortest; ++i) {
        List_Ptr zipper = SASS_MEMORY_NEW(List, pstate, L);
        for (size_t j = 0; j < L; ++j) {
          zipper->append(Cast<List>(arglist->value_at_index(j))->at(i));
        }
        zippers->append(zipper);
      }
      return zippers;
    }

    Signature list_separator_sig = "list_separator($list)";
    BUILT_IN(list_separator)
    {
      List_Obj l = Cast<List>(env["$list"]);
      if (!l) {
        l = SASS_MEMORY_NEW(List, pstate, 1);
        l->append(ARG("$list", Expression));
      }
      return SASS_MEMORY_NEW(String_Quoted,
                               pstate,
                               l->separator() == SASS_COMMA ? "comma" : "space");
    }

    Signature is_bracketed_sig = "is-bracketed($list)";
    BUILT_IN(is_bracketed)
    {
      Value_Obj value = ARG("$list", Value);
      List_Obj list = Cast<List>(value);
      return SASS_MEMORY_NEW(Boolean, pstate, list && list->is_bracketed());
    }

  }

}

namespace Sass {

  Emitter::Emitter(struct Sass_Output_Options& opt)
  : wbuf(),
    opt(opt),
    indentation(0),
    scheduled_space(0),
    scheduled_linefeed(0),
    scheduled_delimiter(false),
    scheduled_crutch(0),
    scheduled_mapping(0),
    in_custom_property(false),
    in_comment(false),
    in_wrapped(false),
    in_media_block(false),
    in_declaration(false),
    in_space_array(false),
    in_comma_array(false)
  { }

  // return buffer as string
  std::string Emitter::get_buffer(void)
  {
    return wbuf.buffer;
  }

  Sass_Output_Style Emitter::output_style(void) const
  {
    return opt.output_style;
  }

  // PROXY METHODS FOR SOURCE MAPS

  void Emitter::add_source_index(size_t idx)
  { wbuf.smap.source_index.push_back(idx); }

  std::string Emitter::render_srcmap(Context &ctx)
  { return wbuf.smap.render_srcmap(ctx); }

  void Emitter::set_filename(const std::string& str)
  { wbuf.smap.file = str; }

  void Emitter::schedule_mapping(AST_Node_Ptr_Const node)
  { scheduled_mapping = node; }
  void Emitter::add_open_mapping(AST_Node_Ptr_Const node)
  { wbuf.smap.add_open_mapping(node); }
  void Emitter::add_close_mapping(AST_Node_Ptr_Const node)
  { wbuf.smap.add_close_mapping(node); }
  ParserState Emitter::remap(const ParserState& pstate)
  { return wbuf.smap.remap(pstate); }

  // MAIN BUFFER MANIPULATION

  // add outstanding delimiter
  void Emitter::finalize(bool final)
  {
    scheduled_space = 0;
    if (output_style() == SASS_STYLE_COMPRESSED)
      if (final) scheduled_delimiter = false;
    if (scheduled_linefeed)
      scheduled_linefeed = 1;
    flush_schedules();
  }

  // flush scheduled space/linefeed
  void Emitter::flush_schedules(void)
  {
    // check the schedule
    if (scheduled_linefeed) {
      std::string linefeeds = "";

      for (size_t i = 0; i < scheduled_linefeed; i++)
        linefeeds += opt.linefeed;
      scheduled_space = 0;
      scheduled_linefeed = 0;
      append_string(linefeeds);

    } else if (scheduled_space) {
      std::string spaces(scheduled_space, ' ');
      scheduled_space = 0;
      append_string(spaces);
    }
    if (scheduled_delimiter) {
      scheduled_delimiter = false;
      append_string(";");
    }
  }

  // prepend some text or token to the buffer
  void Emitter::prepend_output(const OutputBuffer& output)
  {
    wbuf.smap.prepend(output);
    wbuf.buffer = output.buffer + wbuf.buffer;
  }

  // prepend some text or token to the buffer
  void Emitter::prepend_string(const std::string& text)
  {
    // do not adjust mappings for utf8 bom
    // seems they are not counted in any UA
    if (text.compare("\xEF\xBB\xBF") != 0) {
      wbuf.smap.prepend(Offset(text));
    }
    wbuf.buffer = text + wbuf.buffer;
  }

  char Emitter::last_char()
  {
    return wbuf.buffer.back();
  }

  // append a single char to the buffer
  void Emitter::append_char(const char chr)
  {
    // write space/lf
    flush_schedules();
    // add to buffer
    wbuf.buffer += chr;
    // account for data in source-maps
    wbuf.smap.append(Offset(chr));
  }

  // append some text or token to the buffer
  void Emitter::append_string(const std::string& text)
  {

    // write space/lf
    flush_schedules();

    if (in_comment && output_style() == COMPACT) {
      // unescape comment nodes
      std::string out = comment_to_string(text);
      // add to buffer
      wbuf.buffer += out;
      // account for data in source-maps
      wbuf.smap.append(Offset(out));
    } else {
      // add to buffer
      wbuf.buffer += text;
      // account for data in source-maps
      wbuf.smap.append(Offset(text));
    }
  }

  // append some white-space only text
  void Emitter::append_wspace(const std::string& text)
  {
    if (text.empty()) return;
    if (peek_linefeed(text.c_str())) {
      scheduled_space = 0;
      append_mandatory_linefeed();
    }
  }

  // append some text or token to the buffer
  // this adds source-mappings for node start and end
  void Emitter::append_token(const std::string& text, const AST_Node_Ptr node)
  {
    flush_schedules();
    add_open_mapping(node);
    // hotfix for browser issues
    // this is pretty ugly indeed
    if (scheduled_crutch) {
      add_open_mapping(scheduled_crutch);
      scheduled_crutch = 0;
    }
    append_string(text);
    add_close_mapping(node);
  }

  // HELPER METHODS

  void Emitter::append_indentation()
  {
    if (output_style() == COMPRESSED) return;
    if (output_style() == COMPACT) return;
    if (in_declaration && in_comma_array) return;
    if (scheduled_linefeed && indentation)
      scheduled_linefeed = 1;
    std::string indent = "";
    for (size_t i = 0; i < indentation; i++)
      indent += opt.indent;
    append_string(indent);
  }

  void Emitter::append_delimiter()
  {
    scheduled_delimiter = true;
    if (output_style() == COMPACT) {
      if (indentation == 0) {
        append_mandatory_linefeed();
      } else {
        append_mandatory_space();
      }
    } else if (output_style() != COMPRESSED) {
      append_optional_linefeed();
    }
  }

  void Emitter::append_comma_separator()
  {
    // scheduled_space = 0;
    append_string(",");
    append_optional_space();
  }

  void Emitter::append_colon_separator()
  {
    scheduled_space = 0;
    append_string(":");
    if (!in_custom_property) append_optional_space();
  }

  void Emitter::append_mandatory_space()
  {
    scheduled_space = 1;
  }

  void Emitter::append_optional_space()
  {
    if ((output_style() != COMPRESSED) && buffer().size()) {
      unsigned char lst = buffer().at(buffer().length() - 1);
      if (!isspace(lst) || scheduled_delimiter) {
        if (last_char() != '(') {
          append_mandatory_space();
        }
      }
    }
  }

  void Emitter::append_special_linefeed()
  {
    if (output_style() == COMPACT) {
      append_mandatory_linefeed();
      for (size_t p = 0; p < indentation; p++)
        append_string(opt.indent);
    }
  }

  void Emitter::append_optional_linefeed()
  {
    if (in_declaration && in_comma_array) return;
    if (output_style() == COMPACT) {
      append_mandatory_space();
    } else {
      append_mandatory_linefeed();
    }
  }

  void Emitter::append_mandatory_linefeed()
  {
    if (output_style() != COMPRESSED) {
      scheduled_linefeed = 1;
      scheduled_space = 0;
      // flush_schedules();
    }
  }

  void Emitter::append_scope_opener(AST_Node_Ptr node)
  {
    scheduled_linefeed = 0;
    append_optional_space();
    flush_schedules();
    if (node) add_open_mapping(node);
    append_string("{");
    append_optional_linefeed();
    // append_optional_space();
    ++ indentation;
  }
  void Emitter::append_scope_closer(AST_Node_Ptr node)
  {
    -- indentation;
    scheduled_linefeed = 0;
    if (output_style() == COMPRESSED)
      scheduled_delimiter = false;
    if (output_style() == EXPANDED) {
      append_optional_linefeed();
      append_indentation();
    } else {
      append_optional_space();
    }
    append_string("}");
    if (node) add_close_mapping(node);
    append_optional_linefeed();
    if (indentation != 0) return;
    if (output_style() != COMPRESSED)
      scheduled_linefeed = 2;
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  Cssize::Cssize(Context& ctx)
  : ctx(ctx),
    traces(ctx.traces),
    block_stack(BlockStack()),
    p_stack(std::vector<Statement_Ptr>())
  { }

  Statement_Ptr Cssize::parent()
  {
    return p_stack.size() ? p_stack.back() : block_stack.front();
  }

  Block_Ptr Cssize::operator()(Block_Ptr b)
  {
    Block_Obj bb = SASS_MEMORY_NEW(Block, b->pstate(), b->length(), b->is_root());
    // bb->tabs(b->tabs());
    block_stack.push_back(bb);
    append_block(b, bb);
    block_stack.pop_back();
    return bb.detach();
  }

  Statement_Ptr Cssize::operator()(Trace_Ptr t)
  {
    traces.push_back(Backtrace(t->pstate()));
    auto result = t->block()->perform(this);
    traces.pop_back();
    return result;
  }

  Statement_Ptr Cssize::operator()(Declaration_Ptr d)
  {
    String_Obj property = Cast<String>(d->property());

    if (Declaration_Ptr dd = Cast<Declaration>(parent())) {
      String_Obj parent_property = Cast<String>(dd->property());
      property = SASS_MEMORY_NEW(String_Constant,
                                 d->property()->pstate(),
                                 parent_property->to_string() + "-" + property->to_string());
      if (!dd->value()) {
        d->tabs(dd->tabs() + 1);
      }
    }

    Declaration_Obj dd = SASS_MEMORY_NEW(Declaration,
                                      d->pstate(),
                                      property,
                                      d->value(),
                                      d->is_important(),
                                      d->is_custom_property());
    dd->is_indented(d->is_indented());
    dd->tabs(d->tabs());

    p_stack.push_back(dd);
    Block_Obj bb = d->block() ? operator()(d->block()) : NULL;
    p_stack.pop_back();

    if (bb && bb->length()) {
      if (dd->value() && !dd->value()->is_invisible()) {
        bb->unshift(dd);
      }
      return bb.detach();
    }
    else if (dd->value() && !dd->value()->is_invisible()) {
      return dd.detach();
    }

    return 0;
  }

  Statement_Ptr Cssize::operator()(Directive_Ptr r)
  {
    if (!r->block() || !r->block()->length()) return r;

    if (parent()->statement_type() == Statement::RULESET)
    {
      return (r->is_keyframes()) ? SASS_MEMORY_NEW(Bubble, r->pstate(), r) : bubble(r);
    }

    p_stack.push_back(r);
    Directive_Obj rr = SASS_MEMORY_NEW(Directive,
                                  r->pstate(),
                                  r->keyword(),
                                  r->selector(),
                                  r->block() ? operator()(r->block()) : 0);
    if (r->value()) rr->value(r->value());
    p_stack.pop_back();

    bool directive_exists = false;
    size_t L = rr->block() ? rr->block()->length() : 0;
    for (size_t i = 0; i < L && !directive_exists; ++i) {
      Statement_Obj s = r->block()->at(i);
      if (s->statement_type() != Statement::BUBBLE) directive_exists = true;
      else {
        Bubble_Obj s_obj = Cast<Bubble>(s);
        s = s_obj->node();
        if (s->statement_type() != Statement::DIRECTIVE) directive_exists = false;
        else directive_exists = (Cast<Directive>(s)->keyword() == rr->keyword());
      }

    }

    Block_Ptr result = SASS_MEMORY_NEW(Block, rr->pstate());
    if (!(directive_exists || rr->is_keyframes()))
    {
      Directive_Ptr empty_node = Cast<Directive>(rr);
      empty_node->block(SASS_MEMORY_NEW(Block, rr->block() ? rr->block()->pstate() : rr->pstate()));
      result->append(empty_node);
    }

    Block_Obj db = rr->block();
    if (db.isNull()) db = SASS_MEMORY_NEW(Block, rr->pstate());
    Block_Obj ss = debubble(db, rr);
    for (size_t i = 0, L = ss->length(); i < L; ++i) {
      result->append(ss->at(i));
    }

    return result;
  }

  Statement_Ptr Cssize::operator()(Keyframe_Rule_Ptr r)
  {
    if (!r->block() || !r->block()->length()) return r;

    Keyframe_Rule_Obj rr = SASS_MEMORY_NEW(Keyframe_Rule,
                                        r->pstate(),
                                        operator()(r->block()));
    if (!r->name().isNull()) rr->name(r->name());

    return debubble(rr->block(), rr);
  }

  Statement_Ptr Cssize::operator()(Ruleset_Ptr r)
  {
    p_stack.push_back(r);
    // this can return a string schema
    // string schema is not a statement!
    // r->block() is already a string schema
    // and that is comming from propset expand
    Block_Ptr bb = operator()(r->block());
    // this should protect us (at least a bit) from our mess
    // fixing this properly is harder that it should be ...
    if (Cast<Statement>(bb) == NULL) {
      error("Illegal nesting: Only properties may be nested beneath properties.", r->block()->pstate(), traces);
    }
    Ruleset_Obj rr = SASS_MEMORY_NEW(Ruleset,
                                  r->pstate(),
                                  r->selector(),
                                  bb);

    rr->is_root(r->is_root());
    // rr->tabs(r->block()->tabs());
    p_stack.pop_back();

    if (!rr->block()) {
      error("Illegal nesting: Only properties may be nested beneath properties.", r->block()->pstate(), traces);
    }

    Block_Obj props = SASS_MEMORY_NEW(Block, rr->block()->pstate());
    Block_Ptr rules = SASS_MEMORY_NEW(Block, rr->block()->pstate());
    for (size_t i = 0, L = rr->block()->length(); i < L; i++)
    {
      Statement_Ptr s = rr->block()->at(i);
      if (bubblable(s)) rules->append(s);
      if (!bubblable(s)) props->append(s);
    }

    if (props->length())
    {
      Block_Obj pb = SASS_MEMORY_NEW(Block, rr->block()->pstate());
      pb->concat(props);
      rr->block(pb);

      for (size_t i = 0, L = rules->length(); i < L; i++)
      {
        Statement_Ptr stm = rules->at(i);
        stm->tabs(stm->tabs() + 1);
      }

      rules->unshift(rr);
    }

    Block_Ptr ptr = rules;
    rules = debubble(rules);
    void* lp = ptr;
    void* rp = rules;
    if (lp != rp) {
      Block_Obj obj = ptr;
    }

    if (!(!rules->length() ||
          !bubblable(rules->last()) ||
          parent()->statement_type() == Statement::RULESET))
    {
      rules->last()->group_end(true);
    }
    return rules;
  }

  Statement_Ptr Cssize::operator()(Null_Ptr m)
  {
    return 0;
  }

  Statement_Ptr Cssize::operator()(Media_Block_Ptr m)
  {
    if (parent()->statement_type() == Statement::RULESET)
    { return bubble(m); }

    if (parent()->statement_type() == Statement::MEDIA)
    { return SASS_MEMORY_NEW(Bubble, m->pstate(), m); }

    p_stack.push_back(m);

    Media_Block_Obj mm = SASS_MEMORY_NEW(Media_Block,
                                      m->pstate(),
                                      m->media_queries(),
                                      operator()(m->block()));
    mm->tabs(m->tabs());

    p_stack.pop_back();

    return debubble(mm->block(), mm);
  }

  Statement_Ptr Cssize::operator()(Supports_Block_Ptr m)
  {
    if (!m->block()->length())
    { return m; }

    if (parent()->statement_type() == Statement::RULESET)
    { return bubble(m); }

    p_stack.push_back(m);

    Supports_Block_Obj mm = SASS_MEMORY_NEW(Supports_Block,
                                       m->pstate(),
                                       m->condition(),
                                       operator()(m->block()));
    mm->tabs(m->tabs());

    p_stack.pop_back();

    return debubble(mm->block(), mm);
  }

  Statement_Ptr Cssize::operator()(At_Root_Block_Ptr m)
  {
    bool tmp = false;
    for (size_t i = 0, L = p_stack.size(); i < L; ++i) {
      Statement_Ptr s = p_stack[i];
      tmp |= m->exclude_node(s);
    }

    if (!tmp && m->block())
    {
      Block_Ptr bb = operator()(m->block());
      for (size_t i = 0, L = bb->length(); i < L; ++i) {
        // (bb->elements())[i]->tabs(m->tabs());
        Statement_Obj stm = bb->at(i);
        if (bubblable(stm)) stm->tabs(stm->tabs() + m->tabs());
      }
      if (bb->length() && bubblable(bb->last())) bb->last()->group_end(m->group_end());
      return bb;
    }

    if (m->exclude_node(parent()))
    {
      return SASS_MEMORY_NEW(Bubble, m->pstate(), m);
    }

    return bubble(m);
  }

  Statement_Ptr Cssize::bubble(Directive_Ptr m)
  {
    Block_Ptr bb = SASS_MEMORY_NEW(Block, this->parent()->pstate());
    Has_Block_Obj new_rule = Cast<Has_Block>(SASS_MEMORY_COPY(this->parent()));
    new_rule->block(bb);
    new_rule->tabs(this->parent()->tabs());
    new_rule->block()->concat(m->block());

    Block_Obj wrapper_block = SASS_MEMORY_NEW(Block, m->block() ? m->block()->pstate() : m->pstate());
    wrapper_block->append(new_rule);
    Directive_Obj mm = SASS_MEMORY_NEW(Directive,
                                  m->pstate(),
                                  m->keyword(),
                                  m->selector(),
                                  wrapper_block);
    if (m->value()) mm->value(m->value());

    Bubble_Ptr bubble = SASS_MEMORY_NEW(Bubble, mm->pstate(), mm);
    return bubble;
  }

  Statement_Ptr Cssize::bubble(At_Root_Block_Ptr m)
  {
    if (!m || !m->block()) return NULL;
    Block_Ptr bb = SASS_MEMORY_NEW(Block, this->parent()->pstate());
    Has_Block_Obj new_rule = Cast<Has_Block>(SASS_MEMORY_COPY(this->parent()));
    Block_Ptr wrapper_block = SASS_MEMORY_NEW(Block, m->block()->pstate());
    if (new_rule) {
      new_rule->block(bb);
      new_rule->tabs(this->parent()->tabs());
      new_rule->block()->concat(m->block());
      wrapper_block->append(new_rule);
    }

    At_Root_Block_Ptr mm = SASS_MEMORY_NEW(At_Root_Block,
                                        m->pstate(),
                                        wrapper_block,
                                        m->expression());
    Bubble_Ptr bubble = SASS_MEMORY_NEW(Bubble, mm->pstate(), mm);
    return bubble;
  }

  Statement_Ptr Cssize::bubble(Supports_Block_Ptr m)
  {
    Ruleset_Obj parent = Cast<Ruleset>(SASS_MEMORY_COPY(this->parent()));

    Block_Ptr bb = SASS_MEMORY_NEW(Block, parent->block()->pstate());
    Ruleset_Ptr new_rule = SASS_MEMORY_NEW(Ruleset,
                                        parent->pstate(),
                                        parent->selector(),
                                        bb);
    new_rule->tabs(parent->tabs());
    new_rule->block()->concat(m->block());

    Block_Ptr wrapper_block = SASS_MEMORY_NEW(Block, m->block()->pstate());
    wrapper_block->append(new_rule);
    Supports_Block_Ptr mm = SASS_MEMORY_NEW(Supports_Block,
                                       m->pstate(),
                                       m->condition(),
                                       wrapper_block);

    mm->tabs(m->tabs());

    Bubble_Ptr bubble = SASS_MEMORY_NEW(Bubble, mm->pstate(), mm);
    return bubble;
  }

  Statement_Ptr Cssize::bubble(Media_Block_Ptr m)
  {
    Ruleset_Obj parent = Cast<Ruleset>(SASS_MEMORY_COPY(this->parent()));

    Block_Ptr bb = SASS_MEMORY_NEW(Block, parent->block()->pstate());
    Ruleset_Ptr new_rule = SASS_MEMORY_NEW(Ruleset,
                                        parent->pstate(),
                                        parent->selector(),
                                        bb);
    new_rule->tabs(parent->tabs());
    new_rule->block()->concat(m->block());

    Block_Ptr wrapper_block = SASS_MEMORY_NEW(Block, m->block()->pstate());
    wrapper_block->append(new_rule);
    Media_Block_Obj mm = SASS_MEMORY_NEW(Media_Block,
                                      m->pstate(),
                                      m->media_queries(),
                                      wrapper_block);

    mm->tabs(m->tabs());

    return SASS_MEMORY_NEW(Bubble, mm->pstate(), mm);
  }

  bool Cssize::bubblable(Statement_Ptr s)
  {
    return Cast<Ruleset>(s) || s->bubbles();
  }

  Block_Ptr Cssize::flatten(Block_Ptr b)
  {
    Block_Ptr result = SASS_MEMORY_NEW(Block, b->pstate(), 0, b->is_root());
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Ptr ss = b->at(i);
      if (Block_Ptr bb = Cast<Block>(ss)) {
        Block_Obj bs = flatten(bb);
        for (size_t j = 0, K = bs->length(); j < K; ++j) {
          result->append(bs->at(j));
        }
      }
      else {
        result->append(ss);
      }
    }
    return result;
  }

  std::vector<std::pair<bool, Block_Obj>> Cssize::slice_by_bubble(Block_Ptr b)
  {
    std::vector<std::pair<bool, Block_Obj>> results;

    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj value = b->at(i);
      bool key = Cast<Bubble>(value) != NULL;

      if (!results.empty() && results.back().first == key)
      {
        Block_Obj wrapper_block = results.back().second;
        wrapper_block->append(value);
      }
      else
      {
        Block_Ptr wrapper_block = SASS_MEMORY_NEW(Block, value->pstate());
        wrapper_block->append(value);
        results.push_back(std::make_pair(key, wrapper_block));
      }
    }
    return results;
  }

  Block_Ptr Cssize::debubble(Block_Ptr children, Statement_Ptr parent)
  {
    Has_Block_Obj previous_parent;
    std::vector<std::pair<bool, Block_Obj>> baz = slice_by_bubble(children);
    Block_Obj result = SASS_MEMORY_NEW(Block, children->pstate());

    for (size_t i = 0, L = baz.size(); i < L; ++i) {
      bool is_bubble = baz[i].first;
      Block_Obj slice = baz[i].second;

      if (!is_bubble) {
        if (!parent) {
          result->append(slice);
        }
        else if (previous_parent) {
          previous_parent->block()->concat(slice);
        }
        else {
          previous_parent = Cast<Has_Block>(SASS_MEMORY_COPY(parent));
          previous_parent->block(slice);
          previous_parent->tabs(parent->tabs());

          result->append(previous_parent);
        }
        continue;
      }

      for (size_t j = 0, K = slice->length(); j < K; ++j)
      {
        Statement_Ptr ss;
        Statement_Obj stm = slice->at(j);
        // this has to go now here (too bad)
        Bubble_Obj node = Cast<Bubble>(stm);
        Media_Block_Ptr m1 = NULL;
        Media_Block_Ptr m2 = NULL;
        if (parent) m1 = Cast<Media_Block>(parent);
        if (node) m2 = Cast<Media_Block>(node->node());
        if (!parent ||
            parent->statement_type() != Statement::MEDIA ||
            node->node()->statement_type() != Statement::MEDIA ||
            (m1 && m2 && *m1->media_queries() == *m2->media_queries())
          )
        {
          ss = node->node();
        }
        else
        {
          List_Obj mq = merge_media_queries(
            Cast<Media_Block>(node->node()),
            Cast<Media_Block>(parent)
          );
          if (!mq->length()) continue;
          if (Media_Block* b = Cast<Media_Block>(node->node())) {
            b->media_queries(mq);
          }
          ss = node->node();
        }

        if (!ss) continue;

        ss->tabs(ss->tabs() + node->tabs());
        ss->group_end(node->group_end());

        Block_Obj bb = SASS_MEMORY_NEW(Block,
                                    children->pstate(),
                                    children->length(),
                                    children->is_root());
        bb->append(ss->perform(this));

        Block_Obj wrapper_block = SASS_MEMORY_NEW(Block,
                                              children->pstate(),
                                              children->length(),
                                              children->is_root());

        Block_Ptr wrapper = flatten(bb);
        wrapper_block->append(wrapper);

        if (wrapper->length()) {
          previous_parent = {};
        }

        if (wrapper_block) {
          result->append(wrapper_block);
        }
      }
    }

    return flatten(result);
  }

  void Cssize::append_block(Block_Ptr b, Block_Ptr cur)
  {
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj ith = b->at(i)->perform(this);
      if (Block_Ptr bb = Cast<Block>(ith)) {
        for (size_t j = 0, K = bb->length(); j < K; ++j) {
          cur->append(bb->at(j));
        }
      }
      else if (ith) {
        cur->append(ith);
      }
    }
  }

  List_Ptr Cssize::merge_media_queries(Media_Block_Ptr m1, Media_Block_Ptr m2)
  {
    List_Ptr qq = SASS_MEMORY_NEW(List,
                               m1->media_queries()->pstate(),
                               m1->media_queries()->length(),
                               SASS_COMMA);

    for (size_t i = 0, L = m1->media_queries()->length(); i < L; i++) {
      for (size_t j = 0, K = m2->media_queries()->length(); j < K; j++) {
        Expression_Obj l1 = m1->media_queries()->at(i);
        Expression_Obj l2 = m2->media_queries()->at(j);
        Media_Query_Ptr mq1 = Cast<Media_Query>(l1);
        Media_Query_Ptr mq2 = Cast<Media_Query>(l2);
        Media_Query_Ptr mq = merge_media_query(mq1, mq2);
        if (mq) qq->append(mq);
      }
    }

    return qq;
  }


  Media_Query_Ptr Cssize::merge_media_query(Media_Query_Ptr mq1, Media_Query_Ptr mq2)
  {

    std::string type;
    std::string mod;

    std::string m1 = std::string(mq1->is_restricted() ? "only" : mq1->is_negated() ? "not" : "");
    std::string t1 = mq1->media_type() ? mq1->media_type()->to_string(ctx.c_options) : "";
    std::string m2 = std::string(mq2->is_restricted() ? "only" : mq2->is_negated() ? "not" : "");
    std::string t2 = mq2->media_type() ? mq2->media_type()->to_string(ctx.c_options) : "";


    if (t1.empty()) t1 = t2;
    if (t2.empty()) t2 = t1;

    if ((m1 == "not") ^ (m2 == "not")) {
      if (t1 == t2) {
        return 0;
      }
      type = m1 == "not" ? t2 : t1;
      mod = m1 == "not" ? m2 : m1;
    }
    else if (m1 == "not" && m2 == "not") {
      if (t1 != t2) {
        return 0;
      }
      type = t1;
      mod = "not";
    }
    else if (t1 != t2) {
      return 0;
    } else {
      type = t1;
      mod = m1.empty() ? m2 : m1;
    }

    Media_Query_Ptr mm = SASS_MEMORY_NEW(Media_Query,
                                         mq1->pstate(),
                                         {},
                                         mq1->length() + mq2->length(),
                                         mod == "not",
                                         mod == "only");

    if (!type.empty()) {
      mm->media_type(SASS_MEMORY_NEW(String_Quoted, mq1->pstate(), type));
    }

    mm->concat(mq2);
    mm->concat(mq1);

    return mm;
  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {
  // using namespace Lexer;
  using namespace Constants;

  namespace Prelexer {


    /*

        def string_re(open, close)
          /#{open}((?:\\.|\#(?!\{)|[^#{close}\\#])*)(#{close}|#\{)/m
        end
      end

      # A hash of regular expressions that are used for tokenizing strings.
      #
      # The key is a `[Symbol, Boolean]` pair.
      # The symbol represents which style of quotation to use,
      # while the boolean represents whether or not the string
      # is following an interpolated segment.
      STRING_REGULAR_EXPRESSIONS = {
        :double => {
          /#{open}((?:\\.|\#(?!\{)|[^#{close}\\#])*)(#{close}|#\{)/m
          false => string_re('"', '"'),
          true => string_re('', '"')
        },
        :single => {
          false => string_re("'", "'"),
          true => string_re('', "'")
        },
        :uri => {
          false => /url\(#{W}(#{URLCHAR}*?)(#{W}\)|#\{)/,
          true => /(#{URLCHAR}*?)(#{W}\)|#\{)/
        },
        # Defined in https://developer.mozilla.org/en/CSS/@-moz-document as a
        # non-standard version of http://www.w3.org/TR/css3-conditional/
        :url_prefix => {
          false => /url-prefix\(#{W}(#{URLCHAR}*?)(#{W}\)|#\{)/,
          true => /(#{URLCHAR}*?)(#{W}\)|#\{)/
        },
        :domain => {
          false => /domain\(#{W}(#{URLCHAR}*?)(#{W}\)|#\{)/,
          true => /(#{URLCHAR}*?)(#{W}\)|#\{)/
        }
      }
    */

    /*
      /#{open}
        (
          \\.
          |
          \# (?!\{)
          |
          [^#{close}\\#]
        )*
        (#{close}|#\{)
      /m
      false => string_re('"', '"'),
      true => string_re('', '"')
    */
    extern const char string_double_negates[] = "\"\\#";
    const char* re_string_double_close(const char* src)
    {
      return sequence <
        // valid chars
        zero_plus <
          alternatives <
            // escaped char
            sequence <
              exactly <'\\'>,
              any_char
            >,
            // non interpolate hash
            sequence <
              exactly <'#'>,
              negate <
                exactly <'{'>
              >
            >,
            // other valid chars
            neg_class_char <
              string_double_negates
            >
          >
        >,
        // quoted string closer
        // or interpolate opening
        alternatives <
          exactly <'"'>,
          lookahead < exactly< hash_lbrace > >
        >
      >(src);
    }

    const char* re_string_double_open(const char* src)
    {
      return sequence <
        // quoted string opener
        exactly <'"'>,
        // valid chars
        zero_plus <
          alternatives <
            // escaped char
            sequence <
              exactly <'\\'>,
              any_char
            >,
            // non interpolate hash
            sequence <
              exactly <'#'>,
              negate <
                exactly <'{'>
              >
            >,
            // other valid chars
            neg_class_char <
              string_double_negates
            >
          >
        >,
        // quoted string closer
        // or interpolate opening
        alternatives <
          exactly <'"'>,
          lookahead < exactly< hash_lbrace > >
        >
      >(src);
    }

    extern const char string_single_negates[] = "'\\#";
    const char* re_string_single_close(const char* src)
    {
      return sequence <
        // valid chars
        zero_plus <
          alternatives <
            // escaped char
            sequence <
              exactly <'\\'>,
              any_char
            >,
            // non interpolate hash
            sequence <
              exactly <'#'>,
              negate <
                exactly <'{'>
              >
            >,
            // other valid chars
            neg_class_char <
              string_single_negates
            >
          >
        >,
        // quoted string closer
        // or interpolate opening
        alternatives <
          exactly <'\''>,
          lookahead < exactly< hash_lbrace > >
        >
      >(src);
    }

    const char* re_string_single_open(const char* src)
    {
      return sequence <
        // quoted string opener
        exactly <'\''>,
        // valid chars
        zero_plus <
          alternatives <
            // escaped char
            sequence <
              exactly <'\\'>,
              any_char
            >,
            // non interpolate hash
            sequence <
              exactly <'#'>,
              negate <
                exactly <'{'>
              >
            >,
            // other valid chars
            neg_class_char <
              string_single_negates
            >
          >
        >,
        // quoted string closer
        // or interpolate opening
        alternatives <
          exactly <'\''>,
          lookahead < exactly< hash_lbrace > >
        >
      >(src);
    }

    /*
      :uri => {
        false => /url\(#{W}(#{URLCHAR}*?)(#{W}\)|#\{)/,
        true => /(#{URLCHAR}*?)(#{W}\)|#\{)/
      },
    */
    const char* re_string_uri_close(const char* src)
    {
      return sequence <
        non_greedy<
          alternatives<
            class_char< real_uri_chars >,
            uri_character,
            NONASCII,
            ESCAPE
          >,
          alternatives<
            sequence < optional < W >, exactly <')'> >,
            lookahead < exactly< hash_lbrace > >
          >
        >,
        optional <
          sequence < optional < W >, exactly <')'> >
        >
      >(src);
    }

    const char* re_string_uri_open(const char* src)
    {
      return sequence <
        exactly <'u'>,
        exactly <'r'>,
        exactly <'l'>,
        exactly <'('>,
        W,
        alternatives<
          quoted_string,
          non_greedy<
            alternatives<
              class_char< real_uri_chars >,
              uri_character,
              NONASCII,
              ESCAPE
            >,
            alternatives<
              sequence < W, exactly <')'> >,
              exactly< hash_lbrace >
            >
          >
        >
      >(src);
    }

    // Match a line comment (/.*?(?=\n|\r\n?|\Z)/.
    const char* line_comment(const char* src)
    {
      return sequence<
               exactly <
                 slash_slash
               >,
               non_greedy<
                 any_char,
                 end_of_line
               >
             >(src);
    }

    // Match a block comment.
    const char* block_comment(const char* src)
    {
      return sequence<
               delimited_by<
                 slash_star,
                 star_slash,
                 false
               >
             >(src);
    }
    /* not use anymore - remove?
    const char* block_comment_prefix(const char* src) {
      return exactly<slash_star>(src);
    }
    // Match either comment.
    const char* comment(const char* src) {
      return line_comment(src);
    }
    */

    // Match zero plus white-space or line_comments
    const char* optional_css_whitespace(const char* src) {
      return zero_plus< alternatives<spaces, line_comment> >(src);
    }
    const char* css_whitespace(const char* src) {
      return one_plus< alternatives<spaces, line_comment> >(src);
    }
    // Match optional_css_whitepace plus block_comments
    const char* optional_css_comments(const char* src) {
      return zero_plus< alternatives<spaces, line_comment, block_comment> >(src);
    }
    const char* css_comments(const char* src) {
      return one_plus< alternatives<spaces, line_comment, block_comment> >(src);
    }

    // Match one backslash escaped char /\\./
    const char* escape_seq(const char* src)
    {
      return sequence<
        exactly<'\\'>,
        alternatives <
          minmax_range<
            1, 3, xdigit
          >,
          any_char
        >,
        optional <
          exactly <' '>
        >
      >(src);
    }

    // Match identifier start
    const char* identifier_alpha(const char* src)
    {
      return alternatives<
               unicode_seq,
               alpha,
               unicode,
               exactly<'-'>,
               exactly<'_'>,
               NONASCII,
               ESCAPE,
               escape_seq
             >(src);
    }

    // Match identifier after start
    const char* identifier_alnum(const char* src)
    {
      return alternatives<
               unicode_seq,
               alnum,
               unicode,
               exactly<'-'>,
               exactly<'_'>,
               NONASCII,
               ESCAPE,
               escape_seq
             >(src);
    }

    // Match CSS identifiers.
    const char* strict_identifier(const char* src)
    {
      return sequence<
               one_plus < strict_identifier_alpha >,
               zero_plus < strict_identifier_alnum >
               // word_boundary not needed
             >(src);
    }

    // Match CSS identifiers.
    const char* identifier(const char* src)
    {
      return sequence<
               zero_plus< exactly<'-'> >,
               one_plus < identifier_alpha >,
               zero_plus < identifier_alnum >
               // word_boundary not needed
             >(src);
    }

    const char* strict_identifier_alpha(const char* src)
    {
      return alternatives <
               alpha,
               unicode,
               escape_seq,
               exactly<'_'>
             >(src);
    }

    const char* strict_identifier_alnum(const char* src)
    {
      return alternatives <
               alnum,
               unicode,
               escape_seq,
               exactly<'_'>
             >(src);
    }

    // Match a single CSS unit
    const char* one_unit(const char* src)
    {
      return sequence <
               optional < exactly <'-'> >,
               strict_identifier_alpha,
               zero_plus < alternatives<
                 strict_identifier_alnum,
                 sequence <
                   one_plus < exactly<'-'> >,
                   strict_identifier_alpha
                 >
               > >
             >(src);
    }

    // Match numerator/denominator CSS units
    const char* multiple_units(const char* src)
    {
      return
        sequence <
          one_unit,
          zero_plus <
            sequence <
              exactly <'*'>,
              one_unit
            >
          >
        >(src);
    }

    // Match complex CSS unit identifiers
    const char* unit_identifier(const char* src)
    {
      return sequence <
        multiple_units,
        optional <
          sequence <
          exactly <'/'>,
          negate < sequence <
            exactly < calc_fn_kwd >,
            exactly < '(' >
          > >,
          multiple_units
        > >
      >(src);
    }

    const char* identifier_alnums(const char* src)
    {
      return one_plus< identifier_alnum >(src);
    }

    // Match number prefix ([\+\-]+)
    const char* number_prefix(const char* src) {
      return alternatives <
        exactly < '+' >,
        sequence <
          exactly < '-' >,
          optional_css_whitespace,
          exactly< '-' >
        >
      >(src);
    }

    // Match interpolant schemas
    const char* identifier_schema(const char* src) {

      return sequence <
               one_plus <
                 sequence <
                   zero_plus <
                     alternatives <
                       sequence <
                         optional <
                           exactly <'$'>
                         >,
                         identifier
                       >,
                       exactly <'-'>
                     >
                   >,
                   interpolant,
                   zero_plus <
                     alternatives <
                       digits,
                       sequence <
                         optional <
                           exactly <'$'>
                         >,
                         identifier
                       >,
                       quoted_string,
                       exactly<'-'>
                     >
                   >
                 >
               >,
               negate <
                 exactly<'%'>
               >
             > (src);
    }

    // interpolants can be recursive/nested
    const char* interpolant(const char* src) {
      return recursive_scopes< exactly<hash_lbrace>, exactly<rbrace> >(src);
    }

    // $re_squote = /'(?:$re_itplnt|\\.|[^'])*'/
    const char* single_quoted_string(const char* src) {
      // match a single quoted string, while skipping interpolants
      return sequence <
        exactly <'\''>,
        zero_plus <
          alternatives <
            // skip escapes
            sequence <
              exactly < '\\' >,
              re_linebreak
            >,
            escape_seq,
            unicode_seq,
            // skip interpolants
            interpolant,
            // skip non delimiters
            any_char_but < '\'' >
          >
        >,
        exactly <'\''>
      >(src);
    }

    // $re_dquote = /"(?:$re_itp|\\.|[^"])*"/
    const char* double_quoted_string(const char* src) {
      // match a single quoted string, while skipping interpolants
      return sequence <
        exactly <'"'>,
        zero_plus <
          alternatives <
            // skip escapes
            sequence <
              exactly < '\\' >,
              re_linebreak
            >,
            escape_seq,
            unicode_seq,
            // skip interpolants
            interpolant,
            // skip non delimiters
            any_char_but < '"' >
          >
        >,
        exactly <'"'>
      >(src);
    }

    // $re_quoted = /(?:$re_squote|$re_dquote)/
    const char* quoted_string(const char* src) {
      // match a quoted string, while skipping interpolants
      return alternatives<
        single_quoted_string,
        double_quoted_string
      >(src);
    }

    const char* sass_value(const char* src) {
      return alternatives <
        quoted_string,
        identifier,
        percentage,
        hex,
        dimension,
        number
      >(src);
    }

    // this is basically `one_plus < sass_value >`
    // takes care to not parse invalid combinations
    const char* value_combinations(const char* src) {
      // `2px-2px` is invalid combo
      bool was_number = false;
      const char* pos;
      while (src) {
        if ((pos = alternatives < quoted_string, identifier, percentage, hex >(src))) {
          was_number = false;
          src = pos;
        } else if (!was_number && !exactly<'+'>(src) && (pos = alternatives < dimension, number >(src))) {
          was_number = true;
          src = pos;
        } else {
          break;
        }
      }
      return src;
    }

    // must be at least one interpolant
    // can be surrounded by sass values
    // make sure to never parse (dim)(dim)
    // since this wrongly consumes `2px-1px`
    // `2px1px` is valid number (unit `px1px`)
    const char* value_schema(const char* src)
    {
      return sequence <
        one_plus <
          sequence <
            optional < value_combinations >,
            interpolant,
            optional < value_combinations >
          >
        >
      >(src);
    }

    // Match CSS '@' keywords.
    const char* at_keyword(const char* src) {
      return sequence<exactly<'@'>, identifier>(src);
    }

    /*
        tok(%r{
          (
            \\.
          |
            (?!url\()
            [^"'/\#!;\{\}] # "
          |
            /(?![\*\/])
          |
            \#(?!\{)
          |
            !(?![a-z]) # TODO: never consume "!" when issue 1126 is fixed.
          )+
        }xi) || tok(COMMENT) || tok(SINGLE_LINE_COMMENT) || interp_string || interp_uri ||
                interpolation(:warn_for_color)
    */
    const char* re_almost_any_value_token(const char* src) {

      return alternatives <
        one_plus <
          alternatives <
            sequence <
              exactly <'\\'>,
              any_char
            >,
            sequence <
              negate <
                uri_prefix
              >,
              neg_class_char <
                almost_any_value_class
              >
            >,
            sequence <
              exactly <'/'>,
              negate <
                alternatives <
                  exactly <'/'>,
                  exactly <'*'>
                >
              >
            >,
            sequence <
              exactly <'\\'>,
              exactly <'#'>,
              negate <
                exactly <'{'>
              >
            >,
            sequence <
              exactly <'!'>,
              negate <
                alpha
              >
            >
          >
        >,
        block_comment,
        line_comment,
        interpolant,
        space,
        sequence <
          exactly<'u'>,
          exactly<'r'>,
          exactly<'l'>,
          exactly<'('>,
          zero_plus <
            alternatives <
              class_char< real_uri_chars >,
              uri_character,
              NONASCII,
              ESCAPE
            >
          >,
          // false => /url\(#{W}(#{URLCHAR}*?)(#{W}\)|#\{)/,
          // true => /(#{URLCHAR}*?)(#{W}\)|#\{)/
          exactly<')'>
        >
      >(src);
    }

    /*
      DIRECTIVES = Set[:mixin, :include, :function, :return, :debug, :warn, :for,
        :each, :while, :if, :else, :extend, :import, :media, :charset, :content,
        :_moz_document, :at_root, :error]
    */
    const char* re_special_directive(const char* src) {
      return alternatives <
        word < mixin_kwd >,
        word < include_kwd >,
        word < function_kwd >,
        word < return_kwd >,
        word < debug_kwd >,
        word < warn_kwd >,
        word < for_kwd >,
        word < each_kwd >,
        word < while_kwd >,
        word < if_kwd >,
        word < else_kwd >,
        word < extend_kwd >,
        word < import_kwd >,
        word < media_kwd >,
        word < charset_kwd >,
        word < content_kwd >,
        // exactly < moz_document_kwd >,
        word < at_root_kwd >,
        word < error_kwd >
      >(src);
    }

    const char* re_prefixed_directive(const char* src) {
      return sequence <
        optional <
          sequence <
            exactly <'-'>,
            one_plus < alnum >,
            exactly <'-'>
          >
        >,
        exactly < supports_kwd >
      >(src);
    }

    const char* re_reference_combinator(const char* src) {
      return sequence <
        optional <
          sequence <
            zero_plus <
              exactly <'-'>
            >,
            identifier,
            exactly <'|'>
          >
        >,
        zero_plus <
          exactly <'-'>
        >,
        identifier
      >(src);
    }

    const char* static_reference_combinator(const char* src) {
      return sequence <
        exactly <'/'>,
        re_reference_combinator,
        exactly <'/'>
      >(src);
    }

    const char* schema_reference_combinator(const char* src) {
      return sequence <
        exactly <'/'>,
        optional <
          sequence <
            css_ip_identifier,
            exactly <'|'>
          >
        >,
        css_ip_identifier,
        exactly <'/'>
      > (src);
    }

    const char* kwd_import(const char* src) {
      return word<import_kwd>(src);
    }

    const char* kwd_at_root(const char* src) {
      return word<at_root_kwd>(src);
    }

    const char* kwd_with_directive(const char* src) {
      return word<with_kwd>(src);
    }

    const char* kwd_without_directive(const char* src) {
      return word<without_kwd>(src);
    }

    const char* kwd_media(const char* src) {
      return word<media_kwd>(src);
    }

    const char* kwd_supports_directive(const char* src) {
      return word<supports_kwd>(src);
    }

    const char* kwd_mixin(const char* src) {
      return word<mixin_kwd>(src);
    }

    const char* kwd_function(const char* src) {
      return word<function_kwd>(src);
    }

    const char* kwd_return_directive(const char* src) {
      return word<return_kwd>(src);
    }

    const char* kwd_include_directive(const char* src) {
      return word<include_kwd>(src);
    }

    const char* kwd_content_directive(const char* src) {
      return word<content_kwd>(src);
    }

    const char* kwd_charset_directive(const char* src) {
      return word<charset_kwd>(src);
    }

    const char* kwd_extend(const char* src) {
      return word<extend_kwd>(src);
    }


    const char* kwd_if_directive(const char* src) {
      return word<if_kwd>(src);
    }

    const char* kwd_else_directive(const char* src) {
      return word<else_kwd>(src);
    }
    const char* elseif_directive(const char* src) {
      return sequence< exactly< else_kwd >,
                                optional_css_comments,
                                word< if_after_else_kwd > >(src);
    }

    const char* kwd_for_directive(const char* src) {
      return word<for_kwd>(src);
    }

    const char* kwd_from(const char* src) {
      return word<from_kwd>(src);
    }

    const char* kwd_to(const char* src) {
      return word<to_kwd>(src);
    }

    const char* kwd_through(const char* src) {
      return word<through_kwd>(src);
    }

    const char* kwd_each_directive(const char* src) {
      return word<each_kwd>(src);
    }

    const char* kwd_in(const char* src) {
      return word<in_kwd>(src);
    }

    const char* kwd_while_directive(const char* src) {
      return word<while_kwd>(src);
    }

    const char* name(const char* src) {
      return one_plus< alternatives< alnum,
                                     exactly<'-'>,
                                     exactly<'_'>,
                                     escape_seq > >(src);
    }

    const char* kwd_warn(const char* src) {
      return word<warn_kwd>(src);
    }

    const char* kwd_err(const char* src) {
      return word<error_kwd>(src);
    }

    const char* kwd_dbg(const char* src) {
      return word<debug_kwd>(src);
    }

    /* not used anymore - remove?
    const char* directive(const char* src) {
      return sequence< exactly<'@'>, identifier >(src);
    } */

    const char* kwd_null(const char* src) {
      return word<null_kwd>(src);
    }

    const char* css_identifier(const char* src) {
      return sequence <
               zero_plus <
                 exactly <'-'>
               >,
               identifier
             >(src);
    }

    const char* css_ip_identifier(const char* src) {
      return sequence <
               zero_plus <
                 exactly <'-'>
               >,
               alternatives <
                 identifier,
                 interpolant
               >
             >(src);
    }

    // Match CSS type selectors
    const char* namespace_prefix(const char* src) {
      return sequence <
               optional <
                 alternatives <
                   exactly <'*'>,
                   css_identifier
                 >
               >,
               exactly <'|'>,
               negate <
                 exactly <'='>
               >
             >(src);
    }

    // Match CSS type selectors
    const char* namespace_schema(const char* src) {
      return sequence <
               optional <
                 alternatives <
                   exactly <'*'>,
                   css_ip_identifier
                 >
               >,
               exactly<'|'>,
               negate <
                 exactly <'='>
               >
             >(src);
    }

    const char* hyphens_and_identifier(const char* src) {
      return sequence< zero_plus< exactly< '-' > >, identifier_alnums >(src);
    }
    const char* hyphens_and_name(const char* src) {
      return sequence< zero_plus< exactly< '-' > >, name >(src);
    }
    const char* universal(const char* src) {
      return sequence< optional<namespace_schema>, exactly<'*'> >(src);
    }
    // Match CSS id names.
    const char* id_name(const char* src) {
      return sequence<exactly<'#'>, identifier_alnums >(src);
    }
    // Match CSS class names.
    const char* class_name(const char* src) {
      return sequence<exactly<'.'>, identifier >(src);
    }
    // Attribute name in an attribute selector.
    const char* attribute_name(const char* src) {
      return alternatives< sequence< optional<namespace_schema>, identifier>,
                           identifier >(src);
    }
    // match placeholder selectors
    const char* placeholder(const char* src) {
      return sequence<exactly<'%'>, identifier_alnums >(src);
    }
    // Match CSS numeric constants.

    const char* op(const char* src) {
      return class_char<op_chars>(src);
    }
    const char* sign(const char* src) {
      return class_char<sign_chars>(src);
    }
    const char* unsigned_number(const char* src) {
      return alternatives<sequence< zero_plus<digits>,
                                    exactly<'.'>,
                                    one_plus<digits> >,
                          digits>(src);
    }
    const char* number(const char* src) {
      return sequence<
          optional<sign>,
          unsigned_number,
          optional<
            sequence<
              exactly<'e'>,
              optional<sign>,
              unsigned_number
            >
          >
        >(src);
    }
    const char* coefficient(const char* src) {
      return alternatives< sequence< optional<sign>, digits >,
                           sign >(src);
    }
    const char* binomial(const char* src) {
      return sequence <
               optional < sign >,
               optional < digits >,
               exactly <'n'>,
               zero_plus < sequence <
                 optional_css_whitespace, sign,
                 optional_css_whitespace, digits
               > >
             >(src);
    }
    const char* percentage(const char* src) {
      return sequence< number, exactly<'%'> >(src);
    }
    const char* ampersand(const char* src) {
      return exactly<'&'>(src);
    }

    /* not used anymore - remove?
    const char* em(const char* src) {
      return sequence< number, exactly<em_kwd> >(src);
    } */
    const char* dimension(const char* src) {
      return sequence<number, unit_identifier >(src);
    }
    const char* hex(const char* src) {
      const char* p = sequence< exactly<'#'>, one_plus<xdigit> >(src);
      ptrdiff_t len = p - src;
      return (len != 4 && len != 7) ? 0 : p;
    }
    const char* hexa(const char* src) {
      const char* p = sequence< exactly<'#'>, one_plus<xdigit> >(src);
      ptrdiff_t len = p - src;
      return (len != 5 && len != 9) ? 0 : p;
    }
    const char* hex0(const char* src) {
      const char* p = sequence< exactly<'0'>, exactly<'x'>, one_plus<xdigit> >(src);
      ptrdiff_t len = p - src;
      return (len != 5 && len != 8) ? 0 : p;
    }

    /* no longer used - remove?
    const char* rgb_prefix(const char* src) {
      return word<rgb_fn_kwd>(src);
    }*/
    // Match CSS uri specifiers.

    const char* uri_prefix(const char* src) {
      return sequence <
        exactly <
          url_kwd
        >,
        zero_plus <
          sequence <
            exactly <'-'>,
            one_plus <
              alpha
            >
          >
        >,
        exactly <'('>
      >(src);
    }

    // TODO: rename the following two functions
    /* no longer used - remove?
    const char* uri(const char* src) {
      return sequence< exactly<url_kwd>,
                       optional<spaces>,
                       quoted_string,
                       optional<spaces>,
                       exactly<')'> >(src);
    }*/
    /* no longer used - remove?
    const char* url_value(const char* src) {
      return sequence< optional< sequence< identifier, exactly<':'> > >, // optional protocol
                       one_plus< sequence< zero_plus< exactly<'/'> >, filename > >, // one or more folders and/or trailing filename
                       optional< exactly<'/'> > >(src);
    }*/
    /* no longer used - remove?
    const char* url_schema(const char* src) {
      return sequence< optional< sequence< identifier, exactly<':'> > >, // optional protocol
                       filename_schema >(src); // optional trailing slash
    }*/
    // Match CSS "!important" keyword.
    const char* kwd_important(const char* src) {
      return sequence< exactly<'!'>,
                       optional_css_whitespace,
                       word<important_kwd> >(src);
    }
    // Match CSS "!optional" keyword.
    const char* kwd_optional(const char* src) {
      return sequence< exactly<'!'>,
      optional_css_whitespace,
      word<optional_kwd> >(src);
    }
    // Match Sass "!default" keyword.
    const char* default_flag(const char* src) {
      return sequence< exactly<'!'>,
                       optional_css_whitespace,
                       word<default_kwd> >(src);
    }
    // Match Sass "!global" keyword.
    const char* global_flag(const char* src) {
      return sequence< exactly<'!'>,
                       optional_css_whitespace,
                       word<global_kwd> >(src);
    }
    // Match CSS pseudo-class/element prefixes.
    const char* pseudo_prefix(const char* src) {
      return sequence< exactly<':'>, optional< exactly<':'> > >(src);
    }
    // Match CSS function call openers.
    const char* functional_schema(const char* src) {
      return sequence <
               one_plus <
                 sequence <
                   zero_plus <
                     alternatives <
                       identifier,
                       exactly <'-'>
                     >
                   >,
                   one_plus <
                     sequence <
                       interpolant,
                       alternatives <
                         digits,
                         identifier,
                         exactly<'+'>,
                         exactly<'-'>
                       >
                     >
                   >
                 >
               >,
               negate <
                 exactly <'%'>
               >,
               lookahead <
                 exactly <'('>
               >
             > (src);
    }

    const char* re_nothing(const char* src) {
      return src;
    }

    const char* re_functional(const char* src) {
      return sequence< identifier, optional < block_comment >, exactly<'('> >(src);
    }
    const char* re_pseudo_selector(const char* src) {
      return sequence< identifier, optional < block_comment >, exactly<'('> >(src);
    }
    // Match the CSS negation pseudo-class.
    const char* pseudo_not(const char* src) {
      return word< pseudo_not_fn_kwd >(src);
    }
    // Match CSS 'odd' and 'even' keywords for functional pseudo-classes.
    const char* even(const char* src) {
      return word<even_kwd>(src);
    }
    const char* odd(const char* src) {
      return word<odd_kwd>(src);
    }
    // Match CSS attribute-matching operators.
    const char* exact_match(const char* src) { return exactly<'='>(src); }
    const char* class_match(const char* src) { return exactly<tilde_equal>(src); }
    const char* dash_match(const char* src) { return exactly<pipe_equal>(src); }
    const char* prefix_match(const char* src) { return exactly<caret_equal>(src); }
    const char* suffix_match(const char* src) { return exactly<dollar_equal>(src); }
    const char* substring_match(const char* src) { return exactly<star_equal>(src); }
    // Match CSS combinators.
    /* not used anymore - remove?
    const char* adjacent_to(const char* src) {
      return sequence< optional_spaces, exactly<'+'> >(src);
    }
    const char* precedes(const char* src) {
      return sequence< optional_spaces, exactly<'~'> >(src);
    }
    const char* parent_of(const char* src) {
      return sequence< optional_spaces, exactly<'>'> >(src);
    }
    const char* ancestor_of(const char* src) {
      return sequence< spaces, negate< exactly<'{'> > >(src);
    }*/

    // Match SCSS variable names.
    const char* variable(const char* src) {
      return sequence<exactly<'$'>, identifier>(src);
    }

    // parse `calc`, `-a-calc` and `--b-c-calc`
    // but do not parse `foocalc` or `foo-calc`
    const char* calc_fn_call(const char* src) {
      return sequence <
        optional < sequence <
          hyphens,
          one_plus < sequence <
            strict_identifier,
            hyphens
          > >
        > >,
        exactly < calc_fn_kwd >,
        word_boundary
      >(src);
    }

    // Match Sass boolean keywords.
    const char* kwd_true(const char* src) {
      return word<true_kwd>(src);
    }
    const char* kwd_false(const char* src) {
      return word<false_kwd>(src);
    }
    const char* kwd_only(const char* src) {
      return keyword < only_kwd >(src);
    }
    const char* kwd_and(const char* src) {
      return keyword < and_kwd >(src);
    }
    const char* kwd_or(const char* src) {
      return keyword < or_kwd >(src);
    }
    const char* kwd_not(const char* src) {
      return keyword < not_kwd >(src);
    }
    const char* kwd_eq(const char* src) {
      return exactly<eq>(src);
    }
    const char* kwd_neq(const char* src) {
      return exactly<neq>(src);
    }
    const char* kwd_gt(const char* src) {
      return exactly<gt>(src);
    }
    const char* kwd_gte(const char* src) {
      return exactly<gte>(src);
    }
    const char* kwd_lt(const char* src) {
      return exactly<lt>(src);
    }
    const char* kwd_lte(const char* src) {
      return exactly<lte>(src);
    }
    const char* kwd_using(const char* src) {
      return keyword<using_kwd>(src);
    }

    // match specific IE syntax
    const char* ie_progid(const char* src) {
      return sequence <
        word<progid_kwd>,
        exactly<':'>,
        alternatives< identifier_schema, identifier >,
        zero_plus< sequence<
          exactly<'.'>,
          alternatives< identifier_schema, identifier >
        > >,
        zero_plus < sequence<
          exactly<'('>,
          optional_css_whitespace,
          optional < sequence<
            alternatives< variable, identifier_schema, identifier >,
            optional_css_whitespace,
            exactly<'='>,
            optional_css_whitespace,
            alternatives< variable, identifier_schema, identifier, quoted_string, number, hex, hexa >,
            zero_plus< sequence<
              optional_css_whitespace,
              exactly<','>,
              optional_css_whitespace,
              sequence<
                alternatives< variable, identifier_schema, identifier >,
                optional_css_whitespace,
                exactly<'='>,
                optional_css_whitespace,
                alternatives< variable, identifier_schema, identifier, quoted_string, number, hex, hexa >
              >
            > >
          > >,
          optional_css_whitespace,
          exactly<')'>
        > >
      >(src);
    }
    const char* ie_expression(const char* src) {
      return sequence < word<expression_kwd>, exactly<'('>, skip_over_scopes< exactly<'('>, exactly<')'> > >(src);
    }
    const char* ie_property(const char* src) {
      return alternatives < ie_expression, ie_progid >(src);
    }

    // const char* ie_args(const char* src) {
    //   return sequence< alternatives< ie_keyword_arg, value_schema, quoted_string, interpolant, number, identifier, delimited_by< '(', ')', true> >,
    //                    zero_plus< sequence< optional_css_whitespace, exactly<','>, optional_css_whitespace, alternatives< ie_keyword_arg, value_schema, quoted_string, interpolant, number, identifier, delimited_by<'(', ')', true> > > > >(src);
    // }

    const char* ie_keyword_arg_property(const char* src) {
      return alternatives <
          variable,
          identifier_schema,
          identifier
        >(src);
    }
    const char* ie_keyword_arg_value(const char* src) {
      return alternatives <
          variable,
          identifier_schema,
          identifier,
          quoted_string,
          number,
          hex,
          hexa,
          sequence <
            exactly < '(' >,
            skip_over_scopes <
              exactly < '(' >,
              exactly < ')' >
            >
          >
        >(src);
    }

    const char* ie_keyword_arg(const char* src) {
      return sequence <
        ie_keyword_arg_property,
        optional_css_whitespace,
        exactly<'='>,
        optional_css_whitespace,
        ie_keyword_arg_value
      >(src);
    }

    // Path matching functions.
    /* not used anymore - remove?
    const char* folder(const char* src) {
      return sequence< zero_plus< any_char_except<'/'> >,
                       exactly<'/'> >(src);
    }
    const char* folders(const char* src) {
      return zero_plus< folder >(src);
    }*/
    /* not used anymore - remove?
    const char* chunk(const char* src) {
      char inside_str = 0;
      const char* p = src;
      size_t depth = 0;
      while (true) {
        if (!*p) {
          return 0;
        }
        else if (!inside_str && (*p == '"' || *p == '\'')) {
          inside_str = *p;
        }
        else if (*p == inside_str && *(p-1) != '\\') {
          inside_str = 0;
        }
        else if (*p == '(' && !inside_str) {
          ++depth;
        }
        else if (*p == ')' && !inside_str) {
          if (depth == 0) return p;
          else            --depth;
        }
        ++p;
      }
      // unreachable
      return 0;
    }
    */

    // follow the CSS spec more closely and see if this helps us scan URLs correctly
    /* not used anymore - remove?
    const char* NL(const char* src) {
      return alternatives< exactly<'\n'>,
                           sequence< exactly<'\r'>, exactly<'\n'> >,
                           exactly<'\r'>,
                           exactly<'\f'> >(src);
    }*/

    const char* H(const char* src) {
      return std::isxdigit(*src) ? src+1 : 0;
    }

    const char* W(const char* src) {
      return zero_plus< alternatives<
        space,
        exactly< '\t' >,
        exactly< '\r' >,
        exactly< '\n' >,
        exactly< '\f' >
      > >(src);
    }

    const char* UUNICODE(const char* src) {
      return sequence< exactly<'\\'>,
                       between<H, 1, 6>,
                       optional< W >
                       >(src);
    }

    const char* NONASCII(const char* src) {
      return nonascii(src);
    }

    const char* ESCAPE(const char* src) {
      return alternatives<
        UUNICODE,
        sequence<
          exactly<'\\'>,
          alternatives<
            NONASCII,
            escapable_character
          >
        >
      >(src);
    }

    const char* list_terminator(const char* src) {
      return alternatives <
        exactly<';'>,
        exactly<'}'>,
        exactly<'{'>,
        exactly<')'>,
        exactly<']'>,
        exactly<':'>,
        end_of_file,
        exactly<ellipsis>,
        default_flag,
        global_flag
      >(src);
    };

    const char* space_list_terminator(const char* src) {
      return alternatives <
        exactly<','>,
        list_terminator
      >(src);
    };


    // const char* real_uri_prefix(const char* src) {
    //   return alternatives<
    //     exactly< url_kwd >,
    //     exactly< url_prefix_kwd >
    //   >(src);
    // }

    const char* real_uri(const char* src) {
      return sequence<
        exactly< url_kwd >,
        exactly< '(' >,
        W,
        real_uri_value,
        exactly< ')' >
      >(src);
    }

    const char* real_uri_suffix(const char* src) {
      return sequence< W, exactly< ')' > >(src);
    }

    const char* real_uri_value(const char* src) {
      return
      sequence<
        non_greedy<
          alternatives<
            class_char< real_uri_chars >,
            uri_character,
            NONASCII,
            ESCAPE
          >,
          alternatives<
            real_uri_suffix,
            exactly< hash_lbrace >
          >
        >
      >
      (src);
    }

    const char* static_string(const char* src) {
      const char* pos = src;
      const char * s = quoted_string(pos);
      Token t(pos, s);
      const unsigned int p = count_interval< interpolant >(t.begin, t.end);
      return (p == 0) ? t.end : 0;
    }

    const char* unicode_seq(const char* src) {
      return sequence <
        alternatives <
          exactly< 'U' >,
          exactly< 'u' >
        >,
        exactly< '+' >,
        padded_token <
          6, xdigit,
          exactly < '?' >
        >
      >(src);
    }

    const char* static_component(const char* src) {
      return alternatives< identifier,
                           static_string,
                           percentage,
                           hex,
                           hexa,
                           exactly<'|'>,
                           // exactly<'+'>,
                           sequence < number, unit_identifier >,
                           number,
                           sequence< exactly<'!'>, word<important_kwd> >
                          >(src);
    }

    const char* static_property(const char* src) {
      return
        sequence <
          zero_plus<
            sequence <
              optional_css_comments,
              alternatives <
                exactly<','>,
                exactly<'('>,
                exactly<')'>,
                kwd_optional,
                quoted_string,
                interpolant,
                identifier,
                percentage,
                dimension,
                variable,
                alnum,
                sequence <
                  exactly <'\\'>,
                  any_char
                >
              >
            >
          >,
          lookahead <
            sequence <
              optional_css_comments,
              alternatives <
                exactly <';'>,
                exactly <'}'>,
                end_of_file
              >
            >
          >
        >(src);
    }

    const char* static_value(const char* src) {
      return sequence< sequence<
                         static_component,
                         zero_plus< identifier >
                       >,
                       zero_plus < sequence<
                                     alternatives<
                                       sequence< optional_spaces, alternatives<
                                         exactly < '/' >,
                                         exactly < ',' >,
                                         exactly < ' ' >
                                       >, optional_spaces >,
                                       spaces
                                     >,
                                     static_component
                       > >,
                       zero_plus < spaces >,
                       alternatives< exactly<';'>, exactly<'}'> >
                      >(src);
    }

    extern const char css_variable_url_negates[] = "()[]{}\"'#/";
    const char* css_variable_value(const char* src) {
      return sequence<
        alternatives<
          sequence<
            negate< exactly< url_fn_kwd > >,
            one_plus< neg_class_char< css_variable_url_negates > >
          >,
          sequence< exactly<'#'>, negate< exactly<'{'> > >,
          sequence< exactly<'/'>, negate< exactly<'*'> > >,
          static_string,
          real_uri,
          block_comment
        >
      >(src);
    }

    extern const char css_variable_url_top_level_negates[] = "()[]{}\"'#/;";
    const char* css_variable_top_level_value(const char* src) {
      return sequence<
        alternatives<
          sequence<
            negate< exactly< url_fn_kwd > >,
            one_plus< neg_class_char< css_variable_url_top_level_negates > >
          >,
          sequence< exactly<'#'>, negate< exactly<'{'> > >,
          sequence< exactly<'/'>, negate< exactly<'*'> > >,
          static_string,
          real_uri,
          block_comment
        >
      >(src);
    }

    const char* parenthese_scope(const char* src) {
      return sequence <
        exactly < '(' >,
        skip_over_scopes <
          exactly < '(' >,
          exactly < ')' >
        >
      >(src);
    }

    const char* re_selector_list(const char* src) {
      return alternatives <
        // partial bem selector
        sequence <
          ampersand,
          one_plus <
            exactly < '-' >
          >,
          word_boundary,
          optional_spaces
        >,
        // main selector matching
        one_plus <
          alternatives <
            // consume whitespace and comments
            spaces, block_comment, line_comment,
            // match `/deep/` selector (pass-trough)
            // there is no functionality for it yet
            schema_reference_combinator,
            // match selector ops /[*&%,\[\]]/
            class_char < selector_lookahead_ops >,
            // match selector combinators /[>+~]/
            class_char < selector_combinator_ops >,
            // match pseudo selectors
            sequence <
              exactly <'('>,
              optional_spaces,
              optional <re_selector_list>,
              optional_spaces,
              exactly <')'>
            >,
            // match attribute compare operators
            alternatives <
              exact_match, class_match, dash_match,
              prefix_match, suffix_match, substring_match
            >,
            // main selector match
            sequence <
              // allow namespace prefix
              optional < namespace_schema >,
              // modifiers prefixes
              alternatives <
                sequence <
                  exactly <'#'>,
                  // not for interpolation
                  negate < exactly <'{'> >
                >,
                // class match
                exactly <'.'>,
                // single or double colon
                sequence <
                  optional < pseudo_prefix >,
                  // fix libsass issue 2376
                  negate < uri_prefix >
                >
              >,
              // accept hypens in token
              one_plus < sequence <
                // can start with hyphens
                zero_plus <
                  sequence <
                    exactly <'-'>,
                    optional_spaces
                  >
                >,
                // now the main token
                alternatives <
                  kwd_optional,
                  exactly <'*'>,
                  quoted_string,
                  interpolant,
                  identifier,
                  variable,
                  percentage,
                  binomial,
                  dimension,
                  alnum
                >
              > >,
              // can also end with hyphens
              zero_plus < exactly<'-'> >
            >
          >
        >
      >(src);
    }

    const char* type_selector(const char* src) {
      return sequence< optional<namespace_schema>, identifier>(src);
    }
    const char* re_type_selector(const char* src) {
      return alternatives< type_selector, universal, dimension, percentage, number, identifier_alnums >(src);
    }
    const char* re_static_expression(const char* src) {
      return sequence< number, optional_spaces, exactly<'/'>, optional_spaces, number >(src);
    }

    // lexer special_fn: these functions cannot be overloaded
    // (/((-[\w-]+-)?(calc|element)|expression|progid:[a-z\.]*)\(/i)
    const char* re_special_fun(const char* src) {

      // match this first as we test prefix hyphens
      if (const char* calc = calc_fn_call(src)) {
        return calc;
      }

      return sequence <
        optional <
          sequence <
            exactly <'-'>,
            one_plus <
              alternatives <
                alpha,
                exactly <'+'>,
                exactly <'-'>
              >
            >
          >
        >,
        alternatives <
          word < expression_kwd >,
          sequence <
            sequence <
              exactly < progid_kwd >,
              exactly <':'>
            >,
            zero_plus <
              alternatives <
                char_range <'a', 'z'>,
                exactly <'.'>
              >
            >
          >
        >
      >(src);
    }

  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  Definition_Ptr make_native_function(Signature sig, Native_Function func, Context& ctx)
  {
    Parser sig_parser = Parser::from_c_str(sig, ctx, ctx.traces, ParserState("[built-in function]"));
    sig_parser.lex<Prelexer::identifier>();
    std::string name(Util::normalize_underscores(sig_parser.lexed));
    Parameters_Obj params = sig_parser.parse_parameters();
    return SASS_MEMORY_NEW(Definition,
                          ParserState("[built-in function]"),
                          sig,
                          name,
                          params,
                          func,
                          false);
  }

  Definition_Ptr make_c_function(Sass_Function_Entry c_func, Context& ctx)
  {
    using namespace Prelexer;

    const char* sig = sass_function_get_signature(c_func);
    Parser sig_parser = Parser::from_c_str(sig, ctx, ctx.traces, ParserState("[c function]"));
    // allow to overload generic callback plus @warn, @error and @debug with custom functions
    sig_parser.lex < alternatives < identifier, exactly <'*'>,
                                    exactly < Constants::warn_kwd >,
                                    exactly < Constants::error_kwd >,
                                    exactly < Constants::debug_kwd >
                  >              >();
    std::string name(Util::normalize_underscores(sig_parser.lexed));
    Parameters_Obj params = sig_parser.parse_parameters();
    return SASS_MEMORY_NEW(Definition,
                          ParserState("[c function]"),
                          sig,
                          name,
                          params,
                          c_func);
  }

  namespace Functions {

    std::string function_name(Signature sig)
    {
      std::string str(sig);
      return str.substr(0, str.find('('));
    }

    Map_Ptr get_arg_m(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces)
    {
      AST_Node_Ptr value = env[argname];
      if (Map_Ptr map = Cast<Map>(value)) return map;
      List_Ptr list = Cast<List>(value);
      if (list && list->length() == 0) {
        return SASS_MEMORY_NEW(Map, pstate, 0);
      }
      return get_arg<Map>(argname, env, sig, pstate, traces);
    }

    double get_arg_r(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces, double lo, double hi)
    {
      Number_Ptr val = get_arg<Number>(argname, env, sig, pstate, traces);
      Number tmpnr(val);
      tmpnr.reduce();
      double v = tmpnr.value();
      if (!(lo <= v && v <= hi)) {
        std::stringstream msg;
        msg << "argument `" << argname << "` of `" << sig << "` must be between ";
        msg << lo << " and " << hi;
        error(msg.str(), pstate, traces);
      }
      return v;
    }

    Number_Ptr get_arg_n(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces)
    {
      Number_Ptr val = get_arg<Number>(argname, env, sig, pstate, traces);
      val = SASS_MEMORY_COPY(val);
      val->reduce();
      return val;
    }

    double get_arg_val(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces)
    {
      Number_Ptr val = get_arg<Number>(argname, env, sig, pstate, traces);
      Number tmpnr(val);
      tmpnr.reduce();
      return tmpnr.value();
    }

    double color_num(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces)
    {
      Number_Ptr val = get_arg<Number>(argname, env, sig, pstate, traces);
      Number tmpnr(val);
      tmpnr.reduce();
      if (tmpnr.unit() == "%") {
        return std::min(std::max(tmpnr.value() * 255 / 100.0, 0.0), 255.0);
      } else {
        return std::min(std::max(tmpnr.value(), 0.0), 255.0);
      }
    }

    double alpha_num(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces) {
      Number_Ptr val = get_arg<Number>(argname, env, sig, pstate, traces);
      Number tmpnr(val);
      tmpnr.reduce();
      if (tmpnr.unit() == "%") {
        return std::min(std::max(tmpnr.value(), 0.0), 100.0);
      } else {
        return std::min(std::max(tmpnr.value(), 0.0), 1.0);
      }
    }

    Selector_List_Obj get_arg_sels(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces, Context& ctx) {
      Expression_Obj exp = ARG(argname, Expression);
      if (exp->concrete_type() == Expression::NULL_VAL) {
        std::stringstream msg;
        msg << argname << ": null is not a valid selector: it must be a string,\n";
        msg << "a list of strings, or a list of lists of strings for `" << function_name(sig) << "'";
        error(msg.str(), exp->pstate(), traces);
      }
      if (String_Constant_Ptr str = Cast<String_Constant>(exp)) {
        str->quote_mark(0);
      }
      std::string exp_src = exp->to_string(ctx.c_options);
      return Parser::parse_selector(exp_src.c_str(), ctx, traces, exp->pstate(), pstate.src, /*allow_parent=*/false);
    }

    Compound_Selector_Obj get_arg_sel(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces, Context& ctx) {
      Expression_Obj exp = ARG(argname, Expression);
      if (exp->concrete_type() == Expression::NULL_VAL) {
        std::stringstream msg;
        msg << argname << ": null is not a string for `" << function_name(sig) << "'";
        error(msg.str(), exp->pstate(), traces);
      }
      if (String_Constant_Ptr str = Cast<String_Constant>(exp)) {
        str->quote_mark(0);
      }
      std::string exp_src = exp->to_string(ctx.c_options);
      Selector_List_Obj sel_list = Parser::parse_selector(exp_src.c_str(), ctx, traces, exp->pstate(), pstate.src, /*allow_parent=*/false);
      if (sel_list->length() == 0) return {};
      Complex_Selector_Obj first = sel_list->first();
      if (!first->tail()) return first->head();
      return first->tail()->head();
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  // helper to convert string list to vector
  std::vector<std::string> list2vec(struct string_list* cur)
  {
    std::vector<std::string> list;
    while (cur) {
      list.push_back(cur->string);
      cur = cur->next;
    }
    return list;
  }

}

extern "C" {
  using namespace Sass;

  // Allocate libsass heap memory
  // Don't forget string termination!
  void* ADDCALL sass_alloc_memory(size_t size)
  {
    void* ptr = malloc(size);
    if (ptr == NULL) {
      std::cerr << "Out of memory.\n";
      exit(EXIT_FAILURE);
    }
    return ptr;
  }

  char* ADDCALL sass_copy_c_string(const char* str)
  {
    size_t len = strlen(str) + 1;
    char* cpy = (char*) sass_alloc_memory(len);
    std::memcpy(cpy, str, len);
    return cpy;
  }

  // Deallocate libsass heap memory
  void ADDCALL sass_free_memory(void* ptr)
  {
    if (ptr) free (ptr);
  }

  // caller must free the returned memory
  char* ADDCALL sass_string_quote (const char *str, const char quote_mark)
  {
    std::string quoted = quote(str, quote_mark);
    return sass_copy_c_string(quoted.c_str());
  }

  // caller must free the returned memory
  char* ADDCALL sass_string_unquote (const char *str)
  {
    std::string unquoted = unquote(str);
    return sass_copy_c_string(unquoted.c_str());
  }

  char* ADDCALL sass_compiler_find_include (const char* file, struct Sass_Compiler* compiler)
  {
    // get the last import entry to get current base directory
    Sass_Import_Entry import = sass_compiler_get_last_import(compiler);
    const std::vector<std::string>& incs = compiler->cpp_ctx->include_paths;
    // create the vector with paths to lookup
    std::vector<std::string> paths(1 + incs.size());
    paths.push_back(File::dir_name(import->abs_path));
    paths.insert( paths.end(), incs.begin(), incs.end() );
    // now resolve the file path relative to lookup paths
    std::string resolved(File::find_include(file, paths));
    return sass_copy_c_string(resolved.c_str());
  }

  char* ADDCALL sass_compiler_find_file (const char* file, struct Sass_Compiler* compiler)
  {
    // get the last import entry to get current base directory
    Sass_Import_Entry import = sass_compiler_get_last_import(compiler);
    const std::vector<std::string>& incs = compiler->cpp_ctx->include_paths;
    // create the vector with paths to lookup
    std::vector<std::string> paths(1 + incs.size());
    paths.push_back(File::dir_name(import->abs_path));
    paths.insert( paths.end(), incs.begin(), incs.end() );
    // now resolve the file path relative to lookup paths
    std::string resolved(File::find_file(file, paths));
    return sass_copy_c_string(resolved.c_str());
  }

  // Make sure to free the returned value!
  // Incs array has to be null terminated!
  // this has the original resolve logic for sass include
  char* ADDCALL sass_find_include (const char* file, struct Sass_Options* opt)
  {
    std::vector<std::string> vec(list2vec(opt->include_paths));
    std::string resolved(File::find_include(file, vec));
    return sass_copy_c_string(resolved.c_str());
  }

  // Make sure to free the returned value!
  // Incs array has to be null terminated!
  char* ADDCALL sass_find_file (const char* file, struct Sass_Options* opt)
  {
    std::vector<std::string> vec(list2vec(opt->include_paths));
    std::string resolved(File::find_file(file, vec));
    return sass_copy_c_string(resolved.c_str());
  }

  // Get compiled libsass version
  const char* ADDCALL libsass_version(void)
  {
    return LIBSASS_VERSION;
  }

  // Get compiled libsass version
  const char* ADDCALL libsass_language_version(void)
  {
    return LIBSASS_LANGUAGE_VERSION;
  }

}

namespace Sass {

  // helper to aid dreaded MSVC debug mode
  char* sass_copy_string(std::string str)
  {
    // In MSVC the following can lead to segfault:
    // sass_copy_c_string(stream.str().c_str());
    // Reason is that the string returned by str() is disposed before
    // sass_copy_c_string is invoked. The string is actually a stack
    // object, so indeed nobody is holding on to it. So it seems
    // perfectly fair to release it right away. So the const char*
    // by c_str will point to invalid memory. I'm not sure if this is
    // the behavior for all compiler, but I'm pretty sure we would
    // have gotten more issues reported if that would be the case.
    // Wrapping it in a functions seems the cleanest approach as the
    // function must hold on to the stack variable until it's done.
    return sass_copy_c_string(str.c_str());
  }

}
/**
 * sass2scss
 * Licensed under the MIT License
 * Copyright (c) Marcel Greter
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_DEPRECATE
#endif

// include library
#include <stdio.h>

///*
//
// src comments: comments in sass syntax (staring with //)
// css comments: multiline comments in css syntax (starting with /*)
//
// KEEP_COMMENT: keep src comments in the resulting css code
// STRIP_COMMENT: strip out all comments (either src or css)
// CONVERT_COMMENT: convert all src comments to css comments
//
//*/

// our own header

// add namespace for c++
namespace Sass
{

	// return the actual prettify value from options
	#define PRETTIFY(converter) (converter.options - (converter.options & 248))
	// query the options integer to check if the option is enables
	#define KEEP_COMMENT(converter) ((converter.options & SASS2SCSS_KEEP_COMMENT) == SASS2SCSS_KEEP_COMMENT)
	#define STRIP_COMMENT(converter) ((converter.options & SASS2SCSS_STRIP_COMMENT) == SASS2SCSS_STRIP_COMMENT)
	#define CONVERT_COMMENT(converter) ((converter.options & SASS2SCSS_CONVERT_COMMENT) == SASS2SCSS_CONVERT_COMMENT)

	// some makros to access the indentation stack
	#define INDENT(converter) (converter.indents.top())

	// some makros to query comment parser status
	#define IS_PARSING(converter) (converter.comment == "")
	#define IS_COMMENT(converter) (converter.comment != "")
	#define IS_SRC_COMMENT(converter) (converter.comment == "//" && ! CONVERT_COMMENT(converter))
	#define IS_CSS_COMMENT(converter) (converter.comment == "/*" || (converter.comment == "//" && CONVERT_COMMENT(converter)))

	// pretty printer helper function
	static std::string closer (const converter& converter)
	{
		return PRETTIFY(converter) == 0 ? " }" :
		     PRETTIFY(converter) <= 1 ? " }" :
		       "\n" + INDENT(converter) + "}";
	}

	// pretty printer helper function
	static std::string opener (const converter& converter)
	{
		return PRETTIFY(converter) == 0 ? " { " :
		     PRETTIFY(converter) <= 2 ? " {" :
		       "\n" + INDENT(converter) + "{";
	}

	// check if the given string is a pseudo selector
	// needed to differentiate from sass property syntax
	static bool isPseudoSelector (std::string& sel)
	{

		size_t len = sel.length();
		if (len < 1) return false;
		size_t pos = sel.find_first_not_of("abcdefghijklmnopqrstuvwxyz-ABCDEFGHIJKLMNOPQRSTUVWXYZ", 1);
		if (pos != std::string::npos) sel.erase(pos, std::string::npos);
		size_t i = sel.length();
		while (i -- > 0) { sel.at(i) = tolower(sel.at(i)); }

		// CSS Level 1 - Recommendation
		if (sel == ":link") return true;
		if (sel == ":visited") return true;
		if (sel == ":active") return true;

		// CSS Level 2 (Revision 1) - Recommendation
		if (sel == ":lang") return true;
		if (sel == ":first-child") return true;
		if (sel == ":hover") return true;
		if (sel == ":focus") return true;
		// disabled - also valid properties
		// if (sel == ":left") return true;
		// if (sel == ":right") return true;
		if (sel == ":first") return true;

		// Selectors Level 3 - Recommendation
		if (sel == ":target") return true;
		if (sel == ":root") return true;
		if (sel == ":nth-child") return true;
		if (sel == ":nth-last-of-child") return true;
		if (sel == ":nth-of-type") return true;
		if (sel == ":nth-last-of-type") return true;
		if (sel == ":last-child") return true;
		if (sel == ":first-of-type") return true;
		if (sel == ":last-of-type") return true;
		if (sel == ":only-child") return true;
		if (sel == ":only-of-type") return true;
		if (sel == ":empty") return true;
		if (sel == ":not") return true;

		// CSS Basic User Interface Module Level 3 - Working Draft
		if (sel == ":default") return true;
		if (sel == ":valid") return true;
		if (sel == ":invalid") return true;
		if (sel == ":in-range") return true;
		if (sel == ":out-of-range") return true;
		if (sel == ":required") return true;
		if (sel == ":optional") return true;
		if (sel == ":read-only") return true;
		if (sel == ":read-write") return true;
		if (sel == ":dir") return true;
		if (sel == ":enabled") return true;
		if (sel == ":disabled") return true;
		if (sel == ":checked") return true;
		if (sel == ":indeterminate") return true;
		if (sel == ":nth-last-child") return true;

		// Selectors Level 4 - Working Draft
		if (sel == ":any-link") return true;
		if (sel == ":local-link") return true;
		if (sel == ":scope") return true;
		if (sel == ":active-drop-target") return true;
		if (sel == ":valid-drop-target") return true;
		if (sel == ":invalid-drop-target") return true;
		if (sel == ":current") return true;
		if (sel == ":past") return true;
		if (sel == ":future") return true;
		if (sel == ":placeholder-shown") return true;
		if (sel == ":user-error") return true;
		if (sel == ":blank") return true;
		if (sel == ":nth-match") return true;
		if (sel == ":nth-last-match") return true;
		if (sel == ":nth-column") return true;
		if (sel == ":nth-last-column") return true;
		if (sel == ":matches") return true;

		// Fullscreen API - Living Standard
		if (sel == ":fullscreen") return true;

		// not a pseudo selector
		return false;

	}

	static size_t findFirstCharacter (std::string& sass, size_t pos)
	{
		return sass.find_first_not_of(SASS2SCSS_FIND_WHITESPACE, pos);
	}

	static size_t findLastCharacter (std::string& sass, size_t pos)
	{
		return sass.find_last_not_of(SASS2SCSS_FIND_WHITESPACE, pos);
	}

	static bool isUrl (std::string& sass, size_t pos)
	{
		return sass[pos] == 'u' && sass[pos+1] == 'r' && sass[pos+2] == 'l' && sass[pos+3] == '(';
	}

	// check if there is some char data
	// will ignore everything in comments
	static bool hasCharData (std::string& sass)
	{

		size_t col_pos = 0;

		while (true)
		{

			// try to find some meaningfull char
			col_pos = sass.find_first_not_of(" \t\n\v\f\r", col_pos);

			// there was no meaningfull char found
			if (col_pos == std::string::npos) return false;

			// found a multiline comment opener
			if (sass.substr(col_pos, 2) == "/*")
			{
				// find the multiline comment closer
				col_pos = sass.find("*/", col_pos);
				// maybe we did not find the closer here
				if (col_pos == std::string::npos) return false;
				// skip closer
				col_pos += 2;
			}
			else
			{
				return true;
			}

		}

	}
	// EO hasCharData

	// find src comment opener
	// correctly skips quoted strings
	static size_t findCommentOpener (std::string& sass)
	{

		size_t col_pos = 0;
		bool apoed = false;
		bool quoted = false;
		bool comment = false;
		size_t brackets = 0;

		while (col_pos != std::string::npos)
		{

			// process all interesting chars
			col_pos = sass.find_first_of("\"\'/\\*()", col_pos);

			// assertion for valid result
			if (col_pos != std::string::npos)
			{
				char character = sass.at(col_pos);

				if (character == '(')
				{
					if (!quoted && !apoed) brackets ++;
				}
				else if (character == ')')
				{
					if (!quoted && !apoed) brackets --;
				}
				else if (character == '\"')
				{
					// invert quote bool
					if (!apoed && !comment) quoted = !quoted;
				}
				else if (character == '\'')
				{
					// invert quote bool
					if (!quoted && !comment) apoed = !apoed;
				}
				else if (col_pos > 0 && character == '/')
				{
					if (sass.at(col_pos - 1) == '*')
					{
						comment = false;
					}
					// next needs to be a slash too
					else if (sass.at(col_pos - 1) == '/')
					{
						// only found if not in single or double quote, bracket or comment
						if (!quoted && !apoed && !comment && brackets == 0) return col_pos - 1;
					}
				}
				else if (character == '\\')
				{
					// skip next char if in quote
					if (quoted || apoed) col_pos ++;
				}
				// this might be a comment opener
				else if (col_pos > 0 && character == '*')
				{
					// opening a multiline comment
					if (sass.at(col_pos - 1) == '/')
					{
						// we are now in a comment
						if (!quoted && !apoed) comment = true;
					}
				}

				// skip char
				col_pos ++;

			}

		}
		// EO while

		return col_pos;

	}
	// EO findCommentOpener

	// remove multiline comments from sass string
	// correctly skips quoted strings
	static std::string removeMultilineComment (std::string &sass)
	{

		std::string clean = "";
		size_t col_pos = 0;
		size_t open_pos = 0;
		size_t close_pos = 0;
		bool apoed = false;
		bool quoted = false;
		bool comment = false;

		// process sass til string end
		while (col_pos != std::string::npos)
		{

			// process all interesting chars
			col_pos = sass.find_first_of("\"\'/\\*", col_pos);

			// assertion for valid result
			if (col_pos != std::string::npos)
			{
				char character = sass.at(col_pos);

				// found quoted string delimiter
				if (character == '\"')
				{
					if (!apoed && !comment) quoted = !quoted;
				}
				else if (character == '\'')
				{
					if (!quoted && !comment) apoed = !apoed;
				}
				// found possible comment closer
				else if (character == '/')
				{
					// look back to see if it is actually a closer
					if (comment && col_pos > 0 && sass.at(col_pos - 1) == '*')
					{
						close_pos = col_pos + 1; comment = false;
					}
				}
				else if (character == '\\')
				{
					// skip escaped char
					if (quoted || apoed) col_pos ++;
				}
				// this might be a comment opener
				else if (character == '*')
				{
					// look back to see if it is actually an opener
					if (!quoted && !apoed && col_pos > 0 && sass.at(col_pos - 1) == '/')
					{
						comment = true; open_pos = col_pos - 1;
						clean += sass.substr(close_pos, open_pos - close_pos);
					}
				}

				// skip char
				col_pos ++;

			}

		}
		// EO while

		// add final parts (add half open comment text)
		if (comment) clean += sass.substr(open_pos);
		else clean += sass.substr(close_pos);

		// return string
		return clean;

	}
	// EO removeMultilineComment

	// right trim a given string
	std::string rtrim(const std::string &sass)
	{
		std::string trimmed = sass;
		size_t pos_ws = trimmed.find_last_not_of(" \t\n\v\f\r");
		if (pos_ws != std::string::npos)
		{ trimmed.erase(pos_ws + 1); }
		else { trimmed.clear(); }
		return trimmed;
	}
	// EO rtrim

	// flush whitespace and print additional text, but
	// only print additional chars and buffer whitespace
	std::string flush (std::string& sass, converter& converter)
	{

		// return flushed
		std::string scss = "";

		// print whitespace buffer
		scss += PRETTIFY(converter) > 0 ?
		        converter.whitespace : "";
		// reset whitespace buffer
		converter.whitespace = "";

		// remove possible newlines from string
		size_t pos_right = sass.find_last_not_of("\n\r");
		if (pos_right == std::string::npos) return scss;

		// get the linefeeds from the string
		std::string lfs = sass.substr(pos_right + 1);
		sass = sass.substr(0, pos_right + 1);

		// find some source comment opener
		size_t comment_pos = findCommentOpener(sass);
		// check if there was a source comment
		if (comment_pos != std::string::npos)
		{
			// convert comment (but only outside other coments)
			if (CONVERT_COMMENT(converter) && !IS_COMMENT(converter))
			{
				// convert to multiline comment
				sass.at(comment_pos + 1) = '*';
				// add comment node to the whitespace
				sass += " */";
			}
			// not at line start
			if (comment_pos > 0)
			{
				// also include whitespace before the actual comment opener
				size_t ws_pos = sass.find_last_not_of(SASS2SCSS_FIND_WHITESPACE, comment_pos - 1);
				comment_pos = ws_pos == std::string::npos ? 0 : ws_pos + 1;
			}
			if (!STRIP_COMMENT(converter))
			{
				// add comment node to the whitespace
				converter.whitespace += sass.substr(comment_pos);
			}
			else
			{
				// sass = removeMultilineComments(sass);
			}
			// update the actual sass code
			sass = sass.substr(0, comment_pos);
		}

		// add newline as getline discharged it
		converter.whitespace += lfs + "\n";

		// maybe remove any leading whitespace
		if (PRETTIFY(converter) == 0)
		{
			// remove leading whitespace and update string
			size_t pos_left = sass.find_first_not_of(SASS2SCSS_FIND_WHITESPACE);
			if (pos_left != std::string::npos) sass = sass.substr(pos_left);
		}

		// add flushed data
		scss += sass;

		// return string
		return scss;

	}
	// EO flush

	// process a line of the sass text
	std::string process (std::string& sass, converter& converter)
	{

		// resulting string
		std::string scss = "";

		// strip multi line comments
		if (STRIP_COMMENT(converter))
		{
			sass = removeMultilineComment(sass);
		}

		// right trim input
		sass = rtrim(sass);

		// get postion of first meaningfull character in string
		size_t pos_left = sass.find_first_not_of(SASS2SCSS_FIND_WHITESPACE);

		// special case for final run
		if (converter.end_of_file) pos_left = 0;

		// maybe has only whitespace
		if (pos_left == std::string::npos)
		{
			// just add complete whitespace
			converter.whitespace += sass + "\n";
		}
		// have meaningfull first char
		else
		{

			// extract and store indentation string
			std::string indent = sass.substr(0, pos_left);

			// check if current line starts a comment
			std::string open = sass.substr(pos_left, 2);

			// line has less or same indentation
			// finalize previous open parser context
			if (indent.length() <= INDENT(converter).length())
			{

				// close multilinie comment
				if (IS_CSS_COMMENT(converter))
				{
					// check if comments will be stripped anyway
					if (!STRIP_COMMENT(converter)) scss += " */";
				}
				// close src comment comment
				else if (IS_SRC_COMMENT(converter))
				{
					// add a newline to avoid closer on same line
					// this would put the bracket in the comment node
					// no longer needed since we parse them correctly
					// if (KEEP_COMMENT(converter)) scss += "\n";
				}
				// close css properties
				else if (converter.property)
				{
					// add closer unless in concat mode
					if (!converter.comma)
					{
						// if there was no colon we have a selector
						// looks like there were no inner properties
						if (converter.selector) scss += " {}";
						// add final semicolon
						else if (!converter.semicolon) scss += ";";
					}
				}

				// reset comment state
				converter.comment = "";

			}

			// make sure we close every "higher" block
			while (indent.length() < INDENT(converter).length())
			{
				// pop stacked context
				converter.indents.pop();
				// print close bracket
				if (IS_PARSING(converter))
				{ scss += closer(converter); }
				else { scss += " */"; }
				// reset comment state
				converter.comment = "";
			}

			// reset converter state
			converter.selector = false;

			// looks like some undocumented behavior ...
			// https://github.com/mgreter/sass2scss/issues/29
			if (sass.substr(pos_left, 1) == "\\") {
				converter.selector = true;
				sass[pos_left] = ' ';
			}

			// check if we have sass property syntax
			if (sass.substr(pos_left, 1) == ":" && sass.substr(pos_left, 2) != "::")
			{

				// default to a selector
				// change back if property found
				converter.selector = true;
				// get postion of first whitespace char
				size_t pos_wspace = sass.find_first_of(SASS2SCSS_FIND_WHITESPACE, pos_left);
				// assertion check for valid result
				if (pos_wspace != std::string::npos)
				{
					// get the possible pseudo selector
					std::string pseudo = sass.substr(pos_left, pos_wspace - pos_left);
					// get position of the first real property value char
					// pseudo selectors get this far, but have no actual value
					size_t pos_value =  sass.find_first_not_of(SASS2SCSS_FIND_WHITESPACE, pos_wspace);
					// assertion check for valid result
					if (pos_value != std::string::npos)
					{
						// only process if not (fallowed by a semicolon or is a pseudo selector)
						if (!(sass.at(pos_value) == ':' || isPseudoSelector(pseudo)))
						{
							// create new string by interchanging the colon sign for property and value
							sass = indent + sass.substr(pos_left + 1, pos_wspace - pos_left - 1) + ":" + sass.substr(pos_wspace);
							// try to find a colon in the current line, but only ...
							size_t pos_colon = sass.find_first_not_of(":", pos_left);
							// assertion for valid result
							if (pos_colon != std::string::npos)
							{
								// ... after the first word (skip begining colons)
								pos_colon = sass.find_first_of(":", pos_colon);
								// it is a selector if there was no colon found
								converter.selector = pos_colon == std::string::npos;
							}
						}
					}
				}

				// check if we have a BEM property (one colon and no selector)
				if (sass.substr(pos_left, 1) == ":" && converter.selector == true) {
					size_t pos_wspace = sass.find_first_of(SASS2SCSS_FIND_WHITESPACE, pos_left);
					sass = indent + sass.substr(pos_left + 1, pos_wspace) + ":";
				}

			}

			// terminate some statements immediately
			else if (
				sass.substr(pos_left, 5) == "@warn" ||
				sass.substr(pos_left, 6) == "@debug" ||
				sass.substr(pos_left, 6) == "@error" ||
				sass.substr(pos_left, 6) == "@value" ||
				sass.substr(pos_left, 8) == "@charset" ||
				sass.substr(pos_left, 10) == "@namespace"
			) { sass = indent + sass.substr(pos_left); }
			// replace some specific sass shorthand directives (if not fallowed by a white space character)
			else if (sass.substr(pos_left, 1) == "=")
			{ sass = indent + "@mixin " + sass.substr(pos_left + 1); }
			else if (sass.substr(pos_left, 1) == "+")
			{
				// must be followed by a mixin call (no whitespace afterwards or at ending directly)
				if (sass[pos_left+1] != 0 && sass[pos_left+1] != ' ' && sass[pos_left+1] != '\t') {
					sass = indent + "@include " + sass.substr(pos_left + 1);
				}
			}

			// add quotes for import if needed
			else if (sass.substr(pos_left, 7) == "@import")
			{
				// get positions for the actual import url
				size_t pos_import = sass.find_first_of(SASS2SCSS_FIND_WHITESPACE, pos_left + 7);
				size_t pos = sass.find_first_not_of(SASS2SCSS_FIND_WHITESPACE, pos_import);
				size_t start = pos;
				bool in_dqstr = false;
				bool in_sqstr = false;
				bool is_escaped = false;
				do {
					if (is_escaped) {
						is_escaped = false;
					}
					else if (sass[pos] == '\\') {
						is_escaped = true;
					}
					else if (sass[pos] == '"') {
						if (!in_sqstr) in_dqstr = ! in_dqstr;
					}
					else if (sass[pos] == '\'') {
						if (!in_dqstr) in_sqstr = ! in_sqstr;
					}
					else if (in_dqstr || in_sqstr) {
						// skip over quoted stuff
					}
					else if (sass[pos] == ',' || sass[pos] == 0) {
						if (sass[start] != '"' && sass[start] != '\'' && !isUrl(sass, start)) {
							size_t end = findLastCharacter(sass, pos - 1) + 1;
							sass = sass.replace(end, 0, "\"");
							sass = sass.replace(start, 0, "\"");
							pos += 2;
						}
						start = findFirstCharacter(sass, pos + 1);
					}
				}
				while (sass[pos++] != 0);

			}
			else if (
				sass.substr(pos_left, 7) != "@return" &&
				sass.substr(pos_left, 7) != "@extend" &&
				sass.substr(pos_left, 8) != "@include" &&
				sass.substr(pos_left, 8) != "@content"
			) {

				// probably a selector anyway
				converter.selector = true;
				// try to find first colon in the current line
				size_t pos_colon = sass.find_first_of(":", pos_left);
				// assertion that we have a colon
				if (pos_colon != std::string::npos)
				{
					// it is not a selector if we have a space after a colon
					if (sass[pos_colon+1] == ' ') converter.selector = false;
					if (sass[pos_colon+1] == '	') converter.selector = false;
				}

			}

			// current line has more indentation
			if (indent.length() >= INDENT(converter).length())
			{
				// not in comment mode
				if (IS_PARSING(converter))
				{
					// has meaningfull chars
					if (hasCharData(sass))
					{
						// is probably a property
						// also true for selectors
						converter.property = true;
					}
				}
			}
			// current line has more indentation
			if (indent.length() > INDENT(converter).length())
			{
				// not in comment mode
				if (IS_PARSING(converter))
				{
					// had meaningfull chars
					if (converter.property)
					{
						// print block opener
						scss += opener(converter);
						// push new stack context
						converter.indents.push("");
						// store block indentation
						INDENT(converter) = indent;
					}
				}
				// is and will be a src comment
				else if (!IS_CSS_COMMENT(converter))
				{
					// scss does not allow multiline src comments
					// therefore add forward slashes to all lines
					sass.at(INDENT(converter).length()+0) = '/';
					// there is an edge case here if indentation
					// is minimal (will overwrite the fist char)
					sass.at(INDENT(converter).length()+1) = '/';
					// could code around that, but I dont' think
					// this will ever be the cause for any trouble
				}
			}

			// line is opening a new comment
			if (open == "/*" || open == "//")
			{
				// reset the property state
				converter.property = false;
				// close previous comment
				if (IS_CSS_COMMENT(converter) && open != "")
				{
					if (!STRIP_COMMENT(converter) && !CONVERT_COMMENT(converter)) scss += " */";
				}
				// force single line comments
				// into a correct css comment
				if (CONVERT_COMMENT(converter))
				{
					if (IS_PARSING(converter))
					{ sass.at(pos_left + 1) = '*'; }
				}
				// set comment flag
				converter.comment = open;

			}

			// flush data only under certain conditions
			if (!(
				// strip css and src comments if option is set
				(IS_COMMENT(converter) && STRIP_COMMENT(converter)) ||
				// strip src comment even if strip option is not set
				// but only if the keep src comment option is not set
				(IS_SRC_COMMENT(converter) && ! KEEP_COMMENT(converter))
			))
			{
				// flush data and buffer whitespace
				scss += flush(sass, converter);
			}

			// get postion of last meaningfull char
			size_t pos_right = sass.find_last_not_of(SASS2SCSS_FIND_WHITESPACE);

			// check for invalid result
			if (pos_right != std::string::npos)
			{

				// get the last meaningfull char
				std::string close = sass.substr(pos_right, 1);

				// check if next line should be concatenated (list mode)
				converter.comma = IS_PARSING(converter) && close == ",";
				converter.semicolon = IS_PARSING(converter) && close == ";";

				// check if we have more than
				// one meaningfull char
				if (pos_right > 0)
				{

					// get the last two chars from string
					std::string close = sass.substr(pos_right - 1, 2);
					// update parser status for expicitly closed comment
					if (close == "*/") converter.comment = "";

				}

			}
			// EO have meaningfull chars from end

		}
		// EO have meaningfull chars from start

		// return scss
		return scss;

	}
	// EO process

	// read line with either CR, LF or CR LF format
	// http://stackoverflow.com/a/6089413/1550314
	static std::istream& safeGetline(std::istream& is, std::string& t)
	{
		t.clear();

		// The characters in the stream are read one-by-one using a std::streambuf.
		// That is faster than reading them one-by-one using the std::istream.
		// Code that uses streambuf this way must be guarded by a sentry object.
		// The sentry object performs various tasks,
		// such as thread synchronization and updating the stream state.

		std::istream::sentry se(is, true);
		std::streambuf* sb = is.rdbuf();

		for(;;) {
			int c = sb->sbumpc();
			switch (c) {
				case '\n':
					return is;
				case '\r':
					if(sb->sgetc() == '\n')
						sb->sbumpc();
					return is;
				case EOF:
					// Also handle the case when the last line has no line ending
					if(t.empty())
						is.setstate(std::ios::eofbit);
					return is;
				default:
					t += (char)c;
			}
		}
	}

	// the main converter function for c++
	char* sass2scss (const std::string& sass, const int options)
	{

		// local variables
		std::string line;
		std::string scss = "";
		std::stringstream stream(sass);

		// create converter variable
		converter converter;
		// initialise all options
		converter.comma = false;
		converter.property = false;
		converter.selector = false;
		converter.semicolon = false;
		converter.end_of_file = false;
		converter.comment = "";
		converter.whitespace = "";
		converter.indents.push("");
		converter.options = options;

		// read line by line and process them
		while(safeGetline(stream, line) && !stream.eof())
		{ scss += process(line, converter); }

		// create mutable string
		std::string closer = "";
		// set the end of file flag
		converter.end_of_file = true;
		// process to close all open blocks
		scss += process(closer, converter);

		// allocate new memory on the heap
		// caller has to free it after use
		char * cstr = (char*) malloc (scss.length() + 1);
		// create a copy of the string
		strcpy (cstr, scss.c_str());
		// return pointer
		return &cstr[0];

	}
	// EO sass2scss

}
// EO namespace

// implement for c
extern "C"
{

	char* ADDCALL sass2scss (const char* sass, const int options)
	{
		return Sass::sass2scss(sass, options);
	}

	// Get compiled sass2scss version
	const char* ADDCALL sass2scss_version(void) {
		return SASS2SCSS_VERSION;
	}

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {


  /*
    # This is the equivalent of ruby's Sass::Util.paths.
    #
    # Return an array of all possible paths through the given arrays.
    #
    # @param arrs [NodeCollection<NodeCollection<Node>>]
    # @return [NodeCollection<NodeCollection<Node>>]
    #
    # @example
    #   paths([[1, 2], [3, 4], [5]]) #=>
    #     # [[1, 3, 5],
    #     #  [2, 3, 5],
    #     #  [1, 4, 5],
    #     #  [2, 4, 5]]

   The following is the modified version of the ruby code that was more portable to C++. You
   should be able to drop it into ruby 3.2.19 and get the same results from ruby sass.

    def paths(arrs)
        // I changed the inject and maps to an iterative approach to make it easier to implement in C++
      loopStart = [[]]

      for arr in arrs do
        permutations = []
        for e in arr do
          for path in loopStart do
            permutations.push(path + [e])
          end
        end
        loopStart = permutations
      end
    end
  */
  Node paths(const Node& arrs) {

    Node loopStart = Node::createCollection();
    loopStart.collection()->push_back(Node::createCollection());

    for (NodeDeque::iterator arrsIter = arrs.collection()->begin(), arrsEndIter = arrs.collection()->end();
    	arrsIter != arrsEndIter; ++arrsIter) {

      Node& arr = *arrsIter;

      Node permutations = Node::createCollection();

      for (NodeDeque::iterator arrIter = arr.collection()->begin(), arrIterEnd = arr.collection()->end();
      	arrIter != arrIterEnd; ++arrIter) {

        Node& e = *arrIter;

        for (NodeDeque::iterator loopStartIter = loopStart.collection()->begin(), loopStartIterEnd = loopStart.collection()->end();
          loopStartIter != loopStartIterEnd; ++loopStartIter) {

          Node& path = *loopStartIter;

          Node newPermutation = Node::createCollection();
          newPermutation.got_line_feed = arr.got_line_feed;
          newPermutation.plus(path);
          newPermutation.collection()->push_back(e);

          permutations.collection()->push_back(newPermutation);
        }
      }

      loopStart = permutations;
    }

    return loopStart;
  }


  /*
  This is the equivalent of ruby sass' Sass::Util.flatten and [].flatten.
  Sass::Util.flatten requires the number of levels to flatten, while
  [].flatten doesn't and will flatten the entire array. This function
  supports both.

  # Flattens the first `n` nested arrays. If n == -1, all arrays will be flattened
  #
  # @param arr [NodeCollection] The array to flatten
  # @param n [int] The number of levels to flatten
  # @return [NodeCollection] The flattened array

  The following is the modified version of the ruby code that was more portable to C++. You
  should be able to drop it into ruby 3.2.19 and get the same results from ruby sass.

  def flatten(arr, n = -1)
    if n != -1 and n == 0 then
      return arr
    end

    flattened = []

    for e in arr do
      if e.is_a?(Array) then
        flattened.concat(flatten(e, n - 1))
      else
        flattened << e
      end
    end

    return flattened
  end
  */
  Node flatten(Node& arr, int n) {
    if (n != -1 && n == 0) {
      return arr;
    }

    Node flattened = Node::createCollection();
    if (arr.got_line_feed) flattened.got_line_feed = true;

    for (NodeDeque::iterator iter = arr.collection()->begin(), iterEnd = arr.collection()->end();
    	iter != iterEnd; iter++) {
    	Node& e = *iter;

      // e has the lf set
      if (e.isCollection()) {

      	// e.collection().got_line_feed = e.got_line_feed;
      	Node recurseFlattened = flatten(e, n - 1);

      	if(e.got_line_feed) {
      		 flattened.got_line_feed = e.got_line_feed;
      	  recurseFlattened.got_line_feed = e.got_line_feed;
      	}

      	for(auto i : (*recurseFlattened.collection())) {
          if (recurseFlattened.got_line_feed) {

            i.got_line_feed = true;
          }
          flattened.collection()->push_back(i);
      	}

      } else {
      	flattened.collection()->push_back(e);
      }
    }

    return flattened;
  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


extern "C" {
  using namespace Sass;

  Sass_Function_List ADDCALL sass_make_function_list(size_t length)
  {
    return (Sass_Function_List) calloc(length + 1, sizeof(Sass_Function_Entry));
  }

  Sass_Function_Entry ADDCALL sass_make_function(const char* signature, Sass_Function_Fn function, void* cookie)
  {
    Sass_Function_Entry cb = (Sass_Function_Entry) calloc(1, sizeof(Sass_Function));
    if (cb == 0) return 0;
    cb->signature = sass_copy_c_string(signature);
    cb->function = function;
    cb->cookie = cookie;
    return cb;
  }

  void ADDCALL sass_delete_function(Sass_Function_Entry entry)
  {
    free(entry->signature);
    free(entry);
  }

  // Deallocator for the allocated memory
  void ADDCALL sass_delete_function_list(Sass_Function_List list)
  {
    Sass_Function_List it = list;
    if (list == 0) return;
    while(*list) {
      sass_delete_function(*list);
      ++list;
    }
    free(it);
  }

  // Setters and getters for callbacks on function lists
  Sass_Function_Entry ADDCALL sass_function_get_list_entry(Sass_Function_List list, size_t pos) { return list[pos]; }
  void sass_function_set_list_entry(Sass_Function_List list, size_t pos, Sass_Function_Entry cb) { list[pos] = cb; }

  const char* ADDCALL sass_function_get_signature(Sass_Function_Entry cb) { return cb->signature; }
  Sass_Function_Fn ADDCALL sass_function_get_function(Sass_Function_Entry cb) { return cb->function; }
  void* ADDCALL sass_function_get_cookie(Sass_Function_Entry cb) { return cb->cookie; }

  Sass_Importer_Entry ADDCALL sass_make_importer(Sass_Importer_Fn importer, double priority, void* cookie)
  {
    Sass_Importer_Entry cb = (Sass_Importer_Entry) calloc(1, sizeof(Sass_Importer));
    if (cb == 0) return 0;
    cb->importer = importer;
    cb->priority = priority;
    cb->cookie = cookie;
    return cb;
  }

  Sass_Importer_Fn ADDCALL sass_importer_get_function(Sass_Importer_Entry cb) { return cb->importer; }
  double ADDCALL sass_importer_get_priority (Sass_Importer_Entry cb) { return cb->priority; }
  void* ADDCALL sass_importer_get_cookie(Sass_Importer_Entry cb) { return cb->cookie; }

  // Just in case we have some stray import structs
  void ADDCALL sass_delete_importer (Sass_Importer_Entry cb)
  {
    free(cb);
  }

  // Creator for sass custom importer function list
  Sass_Importer_List ADDCALL sass_make_importer_list(size_t length)
  {
    return (Sass_Importer_List) calloc(length + 1, sizeof(Sass_Importer_Entry));
  }

  // Deallocator for the allocated memory
  void ADDCALL sass_delete_importer_list(Sass_Importer_List list)
  {
    Sass_Importer_List it = list;
    if (list == 0) return;
    while(*list) {
      sass_delete_importer(*list);
      ++list;
    }
    free(it);
  }

  Sass_Importer_Entry ADDCALL sass_importer_get_list_entry(Sass_Importer_List list, size_t idx) { return list[idx]; }
  void ADDCALL sass_importer_set_list_entry(Sass_Importer_List list, size_t idx, Sass_Importer_Entry cb) { list[idx] = cb; }

  // Creator for sass custom importer return argument list
  Sass_Import_List ADDCALL sass_make_import_list(size_t length)
  {
    return (Sass_Import**) calloc(length + 1, sizeof(Sass_Import*));
  }

  // Creator for a single import entry returned by the custom importer inside the list
  // We take ownership of the memory for source and srcmap (freed when context is destroyd)
  Sass_Import_Entry ADDCALL sass_make_import(const char* imp_path, const char* abs_path, char* source, char* srcmap)
  {
    Sass_Import* v = (Sass_Import*) calloc(1, sizeof(Sass_Import));
    if (v == 0) return 0;
    v->imp_path = imp_path ? sass_copy_c_string(imp_path) : 0;
    v->abs_path = abs_path ? sass_copy_c_string(abs_path) : 0;
    v->source = source;
    v->srcmap = srcmap;
    v->error = 0;
    v->line = -1;
    v->column = -1;
    return v;
  }

  // Older style, but somehow still valid - keep around or deprecate?
  Sass_Import_Entry ADDCALL sass_make_import_entry(const char* path, char* source, char* srcmap)
  {
    return sass_make_import(path, path, source, srcmap);
  }

  // Upgrade a normal import entry to throw an error (original path can be re-used by error reporting)
  Sass_Import_Entry ADDCALL sass_import_set_error(Sass_Import_Entry import, const char* error, size_t line, size_t col)
  {
    if (import == 0) return 0;
    if (import->error) free(import->error);
    import->error = error ? sass_copy_c_string(error) : 0;
    import->line = line ? line : -1;
    import->column = col ? col : -1;
    return import;
  }

  // Setters and getters for entries on the import list
  void ADDCALL sass_import_set_list_entry(Sass_Import_List list, size_t idx, Sass_Import_Entry entry) { list[idx] = entry; }
  Sass_Import_Entry ADDCALL sass_import_get_list_entry(Sass_Import_List list, size_t idx) { return list[idx]; }

  // Deallocator for the allocated memory
  void ADDCALL sass_delete_import_list(Sass_Import_List list)
  {
    Sass_Import_List it = list;
    if (list == 0) return;
    while(*list) {
      sass_delete_import(*list);
      ++list;
    }
    free(it);
  }

  // Just in case we have some stray import structs
  void ADDCALL sass_delete_import(Sass_Import_Entry import)
  {
    free(import->imp_path);
    free(import->abs_path);
    free(import->source);
    free(import->srcmap);
    free(import->error);
    free(import);
  }

  // Getter for callee entry
  const char* ADDCALL sass_callee_get_name(Sass_Callee_Entry entry) { return entry->name; }
  const char* ADDCALL sass_callee_get_path(Sass_Callee_Entry entry) { return entry->path; }
  size_t ADDCALL sass_callee_get_line(Sass_Callee_Entry entry) { return entry->line; }
  size_t ADDCALL sass_callee_get_column(Sass_Callee_Entry entry) { return entry->column; }
  enum Sass_Callee_Type ADDCALL sass_callee_get_type(Sass_Callee_Entry entry) { return entry->type; }
  Sass_Env_Frame ADDCALL sass_callee_get_env (Sass_Callee_Entry entry) { return &entry->env; }

  // Getters and Setters for environments (lexical, local and global)
  union Sass_Value* ADDCALL sass_env_get_lexical (Sass_Env_Frame env, const char* name) {
    Expression_Ptr ex = Cast<Expression>((*env->frame)[name]);
    return ex != NULL ? ast_node_to_sass_value(ex) : NULL;
  }
  void ADDCALL sass_env_set_lexical (Sass_Env_Frame env, const char* name, union Sass_Value* val) {
    (*env->frame)[name] = sass_value_to_ast_node(val);
  }
  union Sass_Value* ADDCALL sass_env_get_local (Sass_Env_Frame env, const char* name) {
    Expression_Ptr ex = Cast<Expression>(env->frame->get_local(name));
    return ex != NULL ? ast_node_to_sass_value(ex) : NULL;
  }
  void ADDCALL sass_env_set_local (Sass_Env_Frame env, const char* name, union Sass_Value* val) {
    env->frame->set_local(name, sass_value_to_ast_node(val));
  }
  union Sass_Value* ADDCALL sass_env_get_global (Sass_Env_Frame env, const char* name) {
    Expression_Ptr ex = Cast<Expression>(env->frame->get_global(name));
    return ex != NULL ? ast_node_to_sass_value(ex) : NULL;
  }
  void ADDCALL sass_env_set_global (Sass_Env_Frame env, const char* name, union Sass_Value* val) {
    env->frame->set_global(name, sass_value_to_ast_node(val));
  }

  // Getter for import entry
  const char* ADDCALL sass_import_get_imp_path(Sass_Import_Entry entry) { return entry->imp_path; }
  const char* ADDCALL sass_import_get_abs_path(Sass_Import_Entry entry) { return entry->abs_path; }
  const char* ADDCALL sass_import_get_source(Sass_Import_Entry entry) { return entry->source; }
  const char* ADDCALL sass_import_get_srcmap(Sass_Import_Entry entry) { return entry->srcmap; }

  // Getter for import error entry
  size_t ADDCALL sass_import_get_error_line(Sass_Import_Entry entry) { return entry->line; }
  size_t ADDCALL sass_import_get_error_column(Sass_Import_Entry entry) { return entry->column; }
  const char* ADDCALL sass_import_get_error_message(Sass_Import_Entry entry) { return entry->error; }

  // Explicit functions to take ownership of the memory
  // Resets our own property since we do not know if it is still alive
  char* ADDCALL sass_import_take_source(Sass_Import_Entry entry) { char* ptr = entry->source; entry->source = 0; return ptr; }
  char* ADDCALL sass_import_take_srcmap(Sass_Import_Entry entry) { char* ptr = entry->srcmap; entry->srcmap = 0; return ptr; }

}
/*
  Copyright (C) 2011 Joseph A. Adams (joeyadams3.14159@gmail.com)
  All rights reserved.

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#if defined(_MSC_VER) && _MSC_VER < 1900

#include <stdlib.h>
#include <stdarg.h>

static int c99_vsnprintf(char* str, size_t size, const char* format, va_list ap)
{
    int count = -1;

    if (size != 0)
        count = _vsnprintf_s(str, size, _TRUNCATE, format, ap);
    if (count == -1)
        count = _vscprintf(format, ap);

    return count;
}

int snprintf(char* str, size_t size, const char* format, ...)
{
    int count;
    va_list ap;

    va_start(ap, format);
    count = c99_vsnprintf(str, size, format, ap);
    va_end(ap);

    return count;
}

#endif
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {
  SourceMap::SourceMap() : current_position(0, 0, 0), file("stdin") { }
  SourceMap::SourceMap(const std::string& file) : current_position(0, 0, 0), file(file) { }

  std::string SourceMap::render_srcmap(Context &ctx) {

    const bool include_sources = ctx.c_options.source_map_contents;
    const std::vector<std::string> links = ctx.srcmap_links;
    const std::vector<Resource>& sources(ctx.resources);

    JsonNode* json_srcmap = json_mkobject();

    json_append_member(json_srcmap, "version", json_mknumber(3));

    const char *file_name = file.c_str();
    JsonNode *json_file_name = json_mkstring(file_name);
    json_append_member(json_srcmap, "file", json_file_name);

    // pass-through sourceRoot option
    if (!ctx.source_map_root.empty()) {
      JsonNode* root = json_mkstring(ctx.source_map_root.c_str());
      json_append_member(json_srcmap, "sourceRoot", root);
    }

    JsonNode *json_sources = json_mkarray();
    for (size_t i = 0; i < source_index.size(); ++i) {
      std::string source(links[source_index[i]]);
      if (ctx.c_options.source_map_file_urls) {
        source = File::rel2abs(source);
        // check for windows abs path
        if (source[0] == '/') {
          // ends up with three slashes
          source = "file://" + source;
        } else {
          // needs an additional slash
          source = "file:///" + source;
        }
      }
      const char* source_name = source.c_str();
      JsonNode *json_source_name = json_mkstring(source_name);
      json_append_element(json_sources, json_source_name);
    }
    json_append_member(json_srcmap, "sources", json_sources);

    if (include_sources && source_index.size()) {
      JsonNode *json_contents = json_mkarray();
      for (size_t i = 0; i < source_index.size(); ++i) {
        const Resource& resource(sources[source_index[i]]);
        JsonNode *json_content = json_mkstring(resource.contents);
        json_append_element(json_contents, json_content);
      }
      json_append_member(json_srcmap, "sourcesContent", json_contents);
    }

    JsonNode *json_names = json_mkarray();
    // so far we have no implementation for names
    // no problem as we do not alter any identifiers
    json_append_member(json_srcmap, "names", json_names);

    std::string mappings = serialize_mappings();
    JsonNode *json_mappings = json_mkstring(mappings.c_str());
    json_append_member(json_srcmap, "mappings", json_mappings);

    char *str = json_stringify(json_srcmap, "\t");
    std::string result = std::string(str);
    free(str);
    json_delete(json_srcmap);
    return result;
  }

  std::string SourceMap::serialize_mappings() {
    std::string result = "";

    size_t previous_generated_line = 0;
    size_t previous_generated_column = 0;
    size_t previous_original_line = 0;
    size_t previous_original_column = 0;
    size_t previous_original_file = 0;
    for (size_t i = 0; i < mappings.size(); ++i) {
      const size_t generated_line = mappings[i].generated_position.line;
      const size_t generated_column = mappings[i].generated_position.column;
      const size_t original_line = mappings[i].original_position.line;
      const size_t original_column = mappings[i].original_position.column;
      const size_t original_file = mappings[i].original_position.file;

      if (generated_line != previous_generated_line) {
        previous_generated_column = 0;
        if (generated_line > previous_generated_line) {
          result += std::string(generated_line - previous_generated_line, ';');
          previous_generated_line = generated_line;
        }
      }
      else if (i > 0) {
        result += ",";
      }

      // generated column
      result += base64vlq.encode(static_cast<int>(generated_column) - static_cast<int>(previous_generated_column));
      previous_generated_column = generated_column;
      // file
      result += base64vlq.encode(static_cast<int>(original_file) - static_cast<int>(previous_original_file));
      previous_original_file = original_file;
      // source line
      result += base64vlq.encode(static_cast<int>(original_line) - static_cast<int>(previous_original_line));
      previous_original_line = original_line;
      // source column
      result += base64vlq.encode(static_cast<int>(original_column) - static_cast<int>(previous_original_column));
      previous_original_column = original_column;
    }

    return result;
  }

  void SourceMap::prepend(const OutputBuffer& out)
  {
    Offset size(out.smap.current_position);
    for (Mapping mapping : out.smap.mappings) {
      if (mapping.generated_position.line > size.line) {
        throw(std::runtime_error("prepend sourcemap has illegal line"));
      }
      if (mapping.generated_position.line == size.line) {
        if (mapping.generated_position.column > size.column) {
          throw(std::runtime_error("prepend sourcemap has illegal column"));
        }
      }
    }
    // adjust the buffer offset
    prepend(Offset(out.buffer));
    // now add the new mappings
    VECTOR_UNSHIFT(mappings, out.smap.mappings);
  }

  void SourceMap::append(const OutputBuffer& out)
  {
    append(Offset(out.buffer));
  }

  void SourceMap::prepend(const Offset& offset)
  {
    if (offset.line != 0 || offset.column != 0) {
      for (Mapping& mapping : mappings) {
        // move stuff on the first old line
        if (mapping.generated_position.line == 0) {
          mapping.generated_position.column += offset.column;
        }
        // make place for the new lines
        mapping.generated_position.line += offset.line;
      }
    }
    if (current_position.line == 0) {
      current_position.column += offset.column;
    }
    current_position.line += offset.line;
  }

  void SourceMap::append(const Offset& offset)
  {
    current_position += offset;
  }

  void SourceMap::add_open_mapping(AST_Node_Ptr_Const node)
  {
    mappings.push_back(Mapping(node->pstate(), current_position));
  }

  void SourceMap::add_close_mapping(AST_Node_Ptr_Const node)
  {
    mappings.push_back(Mapping(node->pstate() + node->pstate().offset, current_position));
  }

  ParserState SourceMap::remap(const ParserState& pstate) {
    for (size_t i = 0; i < mappings.size(); ++i) {
      if (
        mappings[i].generated_position.file == pstate.file &&
        mappings[i].generated_position.line == pstate.line &&
        mappings[i].generated_position.column == pstate.column
      ) return ParserState(pstate.path, pstate.src, mappings[i].original_position, pstate.offset);
    }
    return ParserState(pstate.path, pstate.src, Position(-1, -1, -1), Offset(0, 0));

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


extern "C" {
  using namespace Sass;

  // Return the sass tag for a generic sass value
  enum Sass_Tag ADDCALL sass_value_get_tag(const union Sass_Value* v) { return v->unknown.tag; }

  // Check value for specified type
  bool ADDCALL sass_value_is_null(const union Sass_Value* v) { return v->unknown.tag == SASS_NULL; }
  bool ADDCALL sass_value_is_number(const union Sass_Value* v) { return v->unknown.tag == SASS_NUMBER; }
  bool ADDCALL sass_value_is_string(const union Sass_Value* v) { return v->unknown.tag == SASS_STRING; }
  bool ADDCALL sass_value_is_boolean(const union Sass_Value* v) { return v->unknown.tag == SASS_BOOLEAN; }
  bool ADDCALL sass_value_is_color(const union Sass_Value* v) { return v->unknown.tag == SASS_COLOR; }
  bool ADDCALL sass_value_is_list(const union Sass_Value* v) { return v->unknown.tag == SASS_LIST; }
  bool ADDCALL sass_value_is_map(const union Sass_Value* v) { return v->unknown.tag == SASS_MAP; }
  bool ADDCALL sass_value_is_error(const union Sass_Value* v) { return v->unknown.tag == SASS_ERROR; }
  bool ADDCALL sass_value_is_warning(const union Sass_Value* v) { return v->unknown.tag == SASS_WARNING; }

  // Getters and setters for Sass_Number
  double ADDCALL sass_number_get_value(const union Sass_Value* v) { return v->number.value; }
  void ADDCALL sass_number_set_value(union Sass_Value* v, double value) { v->number.value = value; }
  const char* ADDCALL sass_number_get_unit(const union Sass_Value* v) { return v->number.unit; }
  void ADDCALL sass_number_set_unit(union Sass_Value* v, char* unit) { v->number.unit = unit; }

  // Getters and setters for Sass_String
  const char* ADDCALL sass_string_get_value(const union Sass_Value* v) { return v->string.value; }
  void ADDCALL sass_string_set_value(union Sass_Value* v, char* value) { v->string.value = value; }
  bool ADDCALL sass_string_is_quoted(const union Sass_Value* v) { return v->string.quoted; }
  void ADDCALL sass_string_set_quoted(union Sass_Value* v, bool quoted) { v->string.quoted = quoted; }

  // Getters and setters for Sass_Boolean
  bool ADDCALL sass_boolean_get_value(const union Sass_Value* v) { return v->boolean.value; }
  void ADDCALL sass_boolean_set_value(union Sass_Value* v, bool value) { v->boolean.value = value; }

  // Getters and setters for Sass_Color
  double ADDCALL sass_color_get_r(const union Sass_Value* v) { return v->color.r; }
  void ADDCALL sass_color_set_r(union Sass_Value* v, double r) { v->color.r = r; }
  double ADDCALL sass_color_get_g(const union Sass_Value* v) { return v->color.g; }
  void ADDCALL sass_color_set_g(union Sass_Value* v, double g) { v->color.g = g; }
  double ADDCALL sass_color_get_b(const union Sass_Value* v) { return v->color.b; }
  void ADDCALL sass_color_set_b(union Sass_Value* v, double b) { v->color.b = b; }
  double ADDCALL sass_color_get_a(const union Sass_Value* v) { return v->color.a; }
  void ADDCALL sass_color_set_a(union Sass_Value* v, double a) { v->color.a = a; }

  // Getters and setters for Sass_List
  size_t ADDCALL sass_list_get_length(const union Sass_Value* v) { return v->list.length; }
  enum Sass_Separator ADDCALL sass_list_get_separator(const union Sass_Value* v) { return v->list.separator; }
  void ADDCALL sass_list_set_separator(union Sass_Value* v, enum Sass_Separator separator) { v->list.separator = separator; }
  bool ADDCALL sass_list_get_is_bracketed(const union Sass_Value* v) { return v->list.is_bracketed; }
  void ADDCALL sass_list_set_is_bracketed(union Sass_Value* v, bool is_bracketed) { v->list.is_bracketed = is_bracketed; }
  // Getters and setters for Sass_List values
  union Sass_Value* ADDCALL sass_list_get_value(const union Sass_Value* v, size_t i) { return v->list.values[i]; }
  void ADDCALL sass_list_set_value(union Sass_Value* v, size_t i, union Sass_Value* value) { v->list.values[i] = value; }

  // Getters and setters for Sass_Map
  size_t ADDCALL sass_map_get_length(const union Sass_Value* v) { return v->map.length; }
  // Getters and setters for Sass_List keys and values
  union Sass_Value* ADDCALL sass_map_get_key(const union Sass_Value* v, size_t i) { return v->map.pairs[i].key; }
  union Sass_Value* ADDCALL sass_map_get_value(const union Sass_Value* v, size_t i) { return v->map.pairs[i].value; }
  void ADDCALL sass_map_set_key(union Sass_Value* v, size_t i, union Sass_Value* key) { v->map.pairs[i].key = key; }
  void ADDCALL sass_map_set_value(union Sass_Value* v, size_t i, union Sass_Value* val) { v->map.pairs[i].value = val; }

  // Getters and setters for Sass_Error
  char* ADDCALL sass_error_get_message(const union Sass_Value* v) { return v->error.message; };
  void ADDCALL sass_error_set_message(union Sass_Value* v, char* msg) { v->error.message = msg; };

  // Getters and setters for Sass_Warning
  char* ADDCALL sass_warning_get_message(const union Sass_Value* v) { return v->warning.message; };
  void ADDCALL sass_warning_set_message(union Sass_Value* v, char* msg) { v->warning.message = msg; };

  // Creator functions for all value types

  union Sass_Value* ADDCALL sass_make_boolean(bool val)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->boolean.tag = SASS_BOOLEAN;
    v->boolean.value = val;
    return v;
  }

  union Sass_Value* ADDCALL sass_make_number(double val, const char* unit)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->number.tag = SASS_NUMBER;
    v->number.value = val;
    v->number.unit = unit ? sass_copy_c_string(unit) : 0;
    if (v->number.unit == 0) { free(v); return 0; }
    return v;
  }

  union Sass_Value* ADDCALL sass_make_color(double r, double g, double b, double a)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->color.tag = SASS_COLOR;
    v->color.r = r;
    v->color.g = g;
    v->color.b = b;
    v->color.a = a;
    return v;
  }

  union Sass_Value* ADDCALL sass_make_string(const char* val)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->string.quoted = false;
    v->string.tag = SASS_STRING;
    v->string.value = val ? sass_copy_c_string(val) : 0;
    if (v->string.value == 0) { free(v); return 0; }
    return v;
  }

  union Sass_Value* ADDCALL sass_make_qstring(const char* val)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->string.quoted = true;
    v->string.tag = SASS_STRING;
    v->string.value = val ? sass_copy_c_string(val) : 0;
    if (v->string.value == 0) { free(v); return 0; }
    return v;
  }

  union Sass_Value* ADDCALL sass_make_list(size_t len, enum Sass_Separator sep, bool is_bracketed)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->list.tag = SASS_LIST;
    v->list.length = len;
    v->list.separator = sep;
    v->list.is_bracketed = is_bracketed;
    v->list.values = (union Sass_Value**) calloc(len, sizeof(union Sass_Value*));
    if (v->list.values == 0) { free(v); return 0; }
    return v;
  }

  union Sass_Value* ADDCALL sass_make_map(size_t len)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->map.tag = SASS_MAP;
    v->map.length = len;
    v->map.pairs = (struct Sass_MapPair*) calloc(len, sizeof(struct Sass_MapPair));
    if (v->map.pairs == 0) { free(v); return 0; }
    return v;
  }

  union Sass_Value* ADDCALL sass_make_null(void)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->null.tag = SASS_NULL;
    return v;
  }

  union Sass_Value* ADDCALL sass_make_error(const char* msg)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->error.tag = SASS_ERROR;
    v->error.message = msg ? sass_copy_c_string(msg) : 0;
    if (v->error.message == 0) { free(v); return 0; }
    return v;
  }

  union Sass_Value* ADDCALL sass_make_warning(const char* msg)
  {
    union Sass_Value* v = (Sass_Value*) calloc(1, sizeof(Sass_Value));
    if (v == 0) return 0;
    v->warning.tag = SASS_WARNING;
    v->warning.message = msg ? sass_copy_c_string(msg) : 0;
    if (v->warning.message == 0) { free(v); return 0; }
    return v;
  }

  // will free all associated sass values
  void ADDCALL sass_delete_value(union Sass_Value* val) {

    size_t i;
    if (val == 0) return;
    switch(val->unknown.tag) {
        case SASS_NULL: {
        }   break;
        case SASS_BOOLEAN: {
        }   break;
        case SASS_NUMBER: {
                free(val->number.unit);
        }   break;
        case SASS_COLOR: {
        }   break;
        case SASS_STRING: {
                free(val->string.value);
        }   break;
        case SASS_LIST: {
                for (i=0; i<val->list.length; i++) {
                    sass_delete_value(val->list.values[i]);
                }
                free(val->list.values);
        }   break;
        case SASS_MAP: {
                for (i=0; i<val->map.length; i++) {
                    sass_delete_value(val->map.pairs[i].key);
                    sass_delete_value(val->map.pairs[i].value);
                }
                free(val->map.pairs);
        }   break;
        case SASS_ERROR: {
                free(val->error.message);
        }   break;
        case SASS_WARNING: {
                free(val->error.message);
        }   break;
        default: break;
    }

    free(val);

    }

  // Make a deep cloned copy of the given sass value
  union Sass_Value* ADDCALL sass_clone_value (const union Sass_Value* val)
  {

    size_t i;
    if (val == 0) return 0;
    switch(val->unknown.tag) {
        case SASS_NULL: {
                return sass_make_null();
        }
        case SASS_BOOLEAN: {
                return sass_make_boolean(val->boolean.value);
        }
        case SASS_NUMBER: {
                return sass_make_number(val->number.value, val->number.unit);
        }
        case SASS_COLOR: {
                return sass_make_color(val->color.r, val->color.g, val->color.b, val->color.a);
        }
        case SASS_STRING: {
                return sass_string_is_quoted(val) ? sass_make_qstring(val->string.value) : sass_make_string(val->string.value);
        }
        case SASS_LIST: {
                union Sass_Value* list = sass_make_list(val->list.length, val->list.separator, val->list.is_bracketed);
                for (i = 0; i < list->list.length; i++) {
                    list->list.values[i] = sass_clone_value(val->list.values[i]);
                }
                return list;
        }
        case SASS_MAP: {
                union Sass_Value* map = sass_make_map(val->map.length);
                for (i = 0; i < val->map.length; i++) {
                    map->map.pairs[i].key = sass_clone_value(val->map.pairs[i].key);
                    map->map.pairs[i].value = sass_clone_value(val->map.pairs[i].value);
                }
                return map;
        }
        case SASS_ERROR: {
                return sass_make_error(val->error.message);
        }
        case SASS_WARNING: {
                return sass_make_warning(val->warning.message);
        }
        default: break;
    }

    return 0;

  }

  union Sass_Value* ADDCALL sass_value_stringify (const union Sass_Value* v, bool compressed, int precision)
  {
    Value_Obj val = sass_value_to_ast_node(v);
    Sass_Inspect_Options options(compressed ? COMPRESSED : NESTED, precision);
    std::string str(val->to_string(options));
    return sass_make_qstring(str.c_str());
  }

  union Sass_Value* ADDCALL sass_value_op (enum Sass_OP op, const union Sass_Value* a, const union Sass_Value* b)
  {

    Sass::Value_Obj rv;

    try {

      Value_Obj lhs = sass_value_to_ast_node(a);
      Value_Obj rhs = sass_value_to_ast_node(b);
      struct Sass_Inspect_Options options(NESTED, 5);

      // see if it's a relational expression
      switch(op) {
        case Sass_OP::EQ:  return sass_make_boolean(Operators::eq(lhs, rhs));
        case Sass_OP::NEQ: return sass_make_boolean(Operators::neq(lhs, rhs));
        case Sass_OP::GT:  return sass_make_boolean(Operators::gt(lhs, rhs));
        case Sass_OP::GTE: return sass_make_boolean(Operators::gte(lhs, rhs));
        case Sass_OP::LT:  return sass_make_boolean(Operators::lt(lhs, rhs));
        case Sass_OP::LTE: return sass_make_boolean(Operators::lte(lhs, rhs));
        case Sass_OP::AND: return ast_node_to_sass_value(lhs->is_false() ? lhs : rhs);
        case Sass_OP::OR:  return ast_node_to_sass_value(lhs->is_false() ? rhs : lhs);
        default: break;
      }

      if (sass_value_is_number(a) && sass_value_is_number(b)) {
        Number_Ptr_Const l_n = Cast<Number>(lhs);
        Number_Ptr_Const r_n = Cast<Number>(rhs);
        rv = Operators::op_numbers(op, *l_n, *r_n, options, l_n->pstate());
      }
      else if (sass_value_is_number(a) && sass_value_is_color(a)) {
        Number_Ptr_Const l_n = Cast<Number>(lhs);
        // Direct HSLA operations are not supported
        // All color maths will be deprecated anyway
        Color_RGBA_Obj r_c = Cast<Color>(rhs)->toRGBA();
        rv = Operators::op_number_color(op, *l_n, *r_c, options, l_n->pstate());
      }
      else if (sass_value_is_color(a) && sass_value_is_number(b)) {
        // Direct HSLA operations are not supported
        // All color maths will be deprecated anyway
        Color_RGBA_Obj l_c = Cast<Color>(lhs)->toRGBA();
        Number_Ptr_Const r_n = Cast<Number>(rhs);
        rv = Operators::op_color_number(op, *l_c, *r_n, options, l_c->pstate());
      }
      else if (sass_value_is_color(a) && sass_value_is_color(b)) {
        // Direct HSLA operations are not supported
        // All color maths will be deprecated anyway
        Color_RGBA_Obj l_c = Cast<Color>(lhs)->toRGBA();
        Color_RGBA_Obj r_c = Cast<Color>(rhs)->toRGBA();
        rv = Operators::op_colors(op, *l_c, *r_c, options, l_c->pstate());
      }
      else /* convert other stuff to string and apply operation */ {
        Value_Ptr l_v = Cast<Value>(lhs);
        Value_Ptr r_v = Cast<Value>(rhs);
        rv = Operators::op_strings(op, *l_v, *r_v, options, l_v->pstate());
      }

      // ToDo: maybe we should should return null value?
      if (!rv) return sass_make_error("invalid return value");

      // convert result back to ast node
      return ast_node_to_sass_value(rv.ptr());
    }

    // simply pass the error message back to the caller for now
    catch (Exception::InvalidSass& e) { return sass_make_error(e.what()); }
    catch (std::bad_alloc&) { return sass_make_error("memory exhausted"); }
    catch (std::exception& e) { return sass_make_error(e.what()); }
    catch (std::string& e) { return sass_make_error(e.c_str()); }
    catch (const char* e) { return sass_make_error(e); }
    catch (...) { return sass_make_error("unknown"); }
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.

#ifdef _WIN32
# ifdef __MINGW32__
#  ifndef off64_t
#   define off64_t _off64_t    /* Workaround for http://sourceforge.net/p/mingw/bugs/2024/ */
#  endif
# endif
# include <direct.h>
# define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
#else
# include <unistd.h>
#endif
#include <fstream>
#include <sys/stat.h>

#ifdef _WIN32
# include <windows.h>

# ifdef _MSC_VER
# include <codecvt>
inline static std::string wstring_to_string(const std::wstring& wstr)
{
    std::wstring_convert<std::codecvt_utf8<wchar_t>, wchar_t> wchar_converter;
    return wchar_converter.to_bytes(wstr);
}
# else // mingw(/gcc) does not support C++11's codecvt yet.
inline static std::string wstring_to_string(const std::wstring &wstr)
{
    int size_needed = WideCharToMultiByte(CP_UTF8, 0, &wstr[0], (int)wstr.size(), NULL, 0, NULL, NULL);
    std::string strTo(size_needed, 0);
    WideCharToMultiByte(CP_UTF8, 0, &wstr[0], (int)wstr.size(), &strTo[0], size_needed, NULL, NULL);
    return strTo;
}
# endif
#endif

namespace Sass {
  namespace File {

    // return the current directory
    // always with forward slashes
    // always with trailing slash
    std::string get_cwd()
    {
      const size_t wd_len = 4096;
      #ifndef _WIN32
        char wd[wd_len];
        char* pwd = getcwd(wd, wd_len);
        // we should check error for more detailed info (e.g. ENOENT)
        // http://man7.org/linux/man-pages/man2/getcwd.2.html#ERRORS
        if (pwd == NULL) throw Exception::OperationError("cwd gone missing");
        std::string cwd = pwd;
      #else
        wchar_t wd[wd_len];
        wchar_t* pwd = _wgetcwd(wd, wd_len);
        if (pwd == NULL) throw Exception::OperationError("cwd gone missing");
        std::string cwd = wstring_to_string(pwd);
        //convert backslashes to forward slashes
        replace(cwd.begin(), cwd.end(), '\\', '/');
      #endif
      if (cwd[cwd.length() - 1] != '/') cwd += '/';
      return cwd;
    }

    // test if path exists and is a file
    bool file_exists(const std::string& path)
    {
      #ifdef _WIN32
        wchar_t resolved[32768];
        // windows unicode filepaths are encoded in utf16
        std::string abspath(join_paths(get_cwd(), path));
        if (!(abspath[0] == '/' && abspath[1] == '/')) {
          abspath = "//?/" + abspath;
        }
        std::wstring wpath(UTF_8::convert_to_utf16(abspath));
        std::replace(wpath.begin(), wpath.end(), '/', '\\');
        DWORD rv = GetFullPathNameW(wpath.c_str(), 32767, resolved, NULL);
        if (rv > 32767) throw Exception::OperationError("Path is too long");
        if (rv == 0) throw Exception::OperationError("Path could not be resolved");
        DWORD dwAttrib = GetFileAttributesW(resolved);
        return (dwAttrib != INVALID_FILE_ATTRIBUTES &&
               (!(dwAttrib & FILE_ATTRIBUTE_DIRECTORY)));
      #else
        struct stat st_buf;
        return (stat (path.c_str(), &st_buf) == 0) &&
               (!S_ISDIR (st_buf.st_mode));
      #endif
    }

    // return if given path is absolute
    // works with *nix and windows paths
    bool is_absolute_path(const std::string& path)
    {
      #ifdef _WIN32
        if (path.length() >= 2 && isalpha(path[0]) && path[1] == ':') return true;
      #endif
      size_t i = 0;
      // check if we have a protocol
      if (path[i] && Prelexer::is_alpha(path[i])) {
        // skip over all alphanumeric characters
        while (path[i] && Prelexer::is_alnum(path[i])) ++i;
        i = i && path[i] == ':' ? i + 1 : 0;
      }
      return path[i] == '/';
    }

    // helper function to find the last directory seperator
    inline size_t find_last_folder_separator(const std::string& path, size_t limit = std::string::npos)
    {
      size_t pos;
      size_t pos_p = path.find_last_of('/', limit);
      #ifdef _WIN32
        size_t pos_w = path.find_last_of('\\', limit);
      #else
        size_t pos_w = std::string::npos;
      #endif
      if (pos_p != std::string::npos && pos_w != std::string::npos) {
        pos = std::max(pos_p, pos_w);
      }
      else if (pos_p != std::string::npos) {
        pos = pos_p;
      }
      else {
        pos = pos_w;
      }
      return pos;
    }

    // return only the directory part of path
    std::string dir_name(const std::string& path)
    {
      size_t pos = find_last_folder_separator(path);
      if (pos == std::string::npos) return "";
      else return path.substr(0, pos+1);
    }

    // return only the filename part of path
    std::string base_name(const std::string& path)
    {
      size_t pos = find_last_folder_separator(path);
      if (pos == std::string::npos) return path;
      else return path.substr(pos+1);
    }

    // do a logical clean up of the path
    // no physical check on the filesystem
    std::string make_canonical_path (std::string path)
    {

      // declarations
      size_t pos;

      #ifdef _WIN32
        //convert backslashes to forward slashes
        replace(path.begin(), path.end(), '\\', '/');
      #endif

      pos = 0; // remove all self references inside the path string
      while((pos = path.find("/./", pos)) != std::string::npos) path.erase(pos, 2);

      // remove all leading and trailing self references
      while(path.length() > 1 && path.substr(0, 2) == "./") path.erase(0, 2);
      while((pos = path.length()) > 1 && path.substr(pos - 2) == "/.") path.erase(pos - 2);


      size_t proto = 0;
      // check if we have a protocol
      if (path[proto] && Prelexer::is_alpha(path[proto])) {
        // skip over all alphanumeric characters
        while (path[proto] && Prelexer::is_alnum(path[proto++])) {}
        // then skip over the mandatory colon
        if (proto && path[proto] == ':') ++ proto;
      }

      // then skip over start slashes
      while (path[proto++] == '/') {}

      pos = proto; // collapse multiple delimiters into a single one
      while((pos = path.find("//", pos)) != std::string::npos) path.erase(pos, 1);

      return path;

    }

    // join two path segments cleanly together
    // but only if right side is not absolute yet
    std::string join_paths(std::string l, std::string r)
    {

      #ifdef _WIN32
        // convert Windows backslashes to URL forward slashes
        replace(l.begin(), l.end(), '\\', '/');
        replace(r.begin(), r.end(), '\\', '/');
      #endif

      if (l.empty()) return r;
      if (r.empty()) return l;

      if (is_absolute_path(r)) return r;
      if (l[l.length()-1] != '/') l += '/';

      // this does a logical cleanup of the right hand path
      // Note that this does collapse x/../y sections into y.
      // This is by design. If /foo on your system is a symlink
      // to /bar/baz, then /foo/../cd is actually /bar/cd,
      // not /cd as a naive ../ removal would give you.
      // will only work on leading double dot dirs on rhs
      // therefore it is safe if lhs is already resolved cwd
      while ((r.length() > 3) && ((r.substr(0, 3) == "../") || (r.substr(0, 3)) == "..\\")) {
        size_t L = l.length(), pos = find_last_folder_separator(l, L - 2);
        bool is_slash = pos + 2 == L && (l[pos+1] == '/' || l[pos+1] == '\\');
        bool is_self = pos + 3 == L && (l[pos+1] == '.');
        if (!is_self && !is_slash) r = r.substr(3);
        else if (pos == std::string::npos) break;
        l = l.substr(0, pos == std::string::npos ? pos : pos + 1);
      }

      return l + r;
    }

    std::string path_for_console(const std::string& rel_path, const std::string& abs_path, const std::string& orig_path)
    {
      // magic algorith goes here!!

      // if the file is outside this directory show the absolute path
      if (rel_path.substr(0, 3) == "../") {
        return orig_path;
      }
      // this seems to work most of the time
      return abs_path == orig_path ? abs_path : rel_path;
    }

    // create an absolute path by resolving relative paths with cwd
    std::string rel2abs(const std::string& path, const std::string& base, const std::string& cwd)
    {
      return make_canonical_path(join_paths(join_paths(cwd + "/", base + "/"), path));
    }

    // create a path that is relative to the given base directory
    // path and base will first be resolved against cwd to make them absolute
    std::string abs2rel(const std::string& path, const std::string& base, const std::string& cwd)
    {

      std::string abs_path = rel2abs(path, cwd);
      std::string abs_base = rel2abs(base, cwd);

      size_t proto = 0;
      // check if we have a protocol
      if (path[proto] && Prelexer::is_alpha(path[proto])) {
        // skip over all alphanumeric characters
        while (path[proto] && Prelexer::is_alnum(path[proto++])) {}
        // then skip over the mandatory colon
        if (proto && path[proto] == ':') ++ proto;
      }

      // distinguish between windows absolute paths and valid protocols
      // we assume that protocols must at least have two chars to be valid
      if (proto && path[proto++] == '/' && proto > 3) return path;

      #ifdef _WIN32
        // absolute link must have a drive letter, and we know that we
        // can only create relative links if both are on the same drive
        if (abs_base[0] != abs_path[0]) return abs_path;
      #endif

      std::string stripped_uri = "";
      std::string stripped_base = "";

      size_t index = 0;
      size_t minSize = std::min(abs_path.size(), abs_base.size());
      for (size_t i = 0; i < minSize; ++i) {
        #ifdef FS_CASE_SENSITIVE
          if (abs_path[i] != abs_base[i]) break;
        #else
          // compare the charactes in a case insensitive manner
          // windows fs is only case insensitive in ascii ranges
          if (tolower(abs_path[i]) != tolower(abs_base[i])) break;
        #endif
        if (abs_path[i] == '/') index = i + 1;
      }
      for (size_t i = index; i < abs_path.size(); ++i) {
        stripped_uri += abs_path[i];
      }
      for (size_t i = index; i < abs_base.size(); ++i) {
        stripped_base += abs_base[i];
      }

      size_t left = 0;
      size_t directories = 0;
      for (size_t right = 0; right < stripped_base.size(); ++right) {
        if (stripped_base[right] == '/') {
          if (stripped_base.substr(left, 2) != "..") {
            ++directories;
          }
          else if (directories > 1) {
            --directories;
          }
          else {
            directories = 0;
          }
          left = right + 1;
        }
      }

      std::string result = "";
      for (size_t i = 0; i < directories; ++i) {
        result += "../";
      }
      result += stripped_uri;

      return result;
    }

    // Resolution order for ambiguous imports:
    // (1) filename as given
    // (2) underscore + given
    // (3) underscore + given + extension
    // (4) given + extension
    // (5) given + _index.scss
    // (6) given + _index.sass
    std::vector<Include> resolve_includes(const std::string& root, const std::string& file, const std::vector<std::string>& exts)
    {
      std::string filename = join_paths(root, file);
      // split the filename
      std::string base(dir_name(file));
      std::string name(base_name(file));
      std::vector<Include> includes;
      // create full path (maybe relative)
      std::string rel_path(join_paths(base, name));
      std::string abs_path(join_paths(root, rel_path));
      if (file_exists(abs_path)) includes.push_back({{ rel_path, root }, abs_path });
      // next test variation with underscore
      rel_path = join_paths(base, "_" + name);
      abs_path = join_paths(root, rel_path);
      if (file_exists(abs_path)) includes.push_back({{ rel_path, root }, abs_path });
      // next test exts plus underscore
      for(auto ext : exts) {
        rel_path = join_paths(base, "_" + name + ext);
        abs_path = join_paths(root, rel_path);
        if (file_exists(abs_path)) includes.push_back({{ rel_path, root }, abs_path });
      }
      // next test plain name with exts
      for(auto ext : exts) {
        rel_path = join_paths(base, name + ext);
        abs_path = join_paths(root, rel_path);
        if (file_exists(abs_path)) includes.push_back({{ rel_path, root }, abs_path });
      }
      // index files
      if (includes.size() == 0) {
        // ignore directories that look like @import'able filename
        for(auto ext : exts) {
          if (ends_with(name, ext)) return includes;
        }
        // next test underscore index exts
        for(auto ext : exts) {
          rel_path = join_paths(base, join_paths(name, "_index" + ext));
          abs_path = join_paths(root, rel_path);
          if (file_exists(abs_path)) includes.push_back({{ rel_path, root }, abs_path });
        }
        // next test plain index exts
        for(auto ext : exts) {
          rel_path = join_paths(base, join_paths(name, "index" + ext));
          abs_path = join_paths(root, rel_path);
          if (file_exists(abs_path)) includes.push_back({{ rel_path, root }, abs_path });
        }
      }
      // nothing found
      return includes;
    }

    std::vector<std::string> find_files(const std::string& file, const std::vector<std::string> paths)
    {
      std::vector<std::string> includes;
      for (std::string path : paths) {
        std::string abs_path(join_paths(path, file));
        if (file_exists(abs_path)) includes.push_back(abs_path);
      }
      return includes;
    }

    std::vector<std::string> find_files(const std::string& file, struct Sass_Compiler* compiler)
    {
      // get the last import entry to get current base directory
      // struct Sass_Options* options = sass_compiler_get_options(compiler);
      Sass_Import_Entry import = sass_compiler_get_last_import(compiler);
      const std::vector<std::string>& incs = compiler->cpp_ctx->include_paths;
      // create the vector with paths to lookup
      std::vector<std::string> paths(1 + incs.size());
      paths.push_back(dir_name(import->abs_path));
      paths.insert(paths.end(), incs.begin(), incs.end());
      // dispatch to find files in paths
      return find_files(file, paths);
    }

    // helper function to search one file in all include paths
    // this is normally not used internally by libsass (C-API sugar)
    std::string find_file(const std::string& file, const std::vector<std::string> paths)
    {
      if (file.empty()) return file;
      auto res = find_files(file, paths);
      return res.empty() ? "" : res.front();
    }

    // helper function to resolve a filename
    std::string find_include(const std::string& file, const std::vector<std::string> paths)
    {
      // search in every include path for a match
      for (size_t i = 0, S = paths.size(); i < S; ++i)
      {
        std::vector<Include> resolved(resolve_includes(paths[i], file));
        if (resolved.size()) return resolved[0].abs_path;
      }
      // nothing found
      return std::string("");
    }

    // try to load the given filename
    // returned memory must be freed
    // will auto convert .sass files
    char* read_file(const std::string& path)
    {
      #ifdef _WIN32
        BYTE* pBuffer;
        DWORD dwBytes;
        wchar_t resolved[32768];
        // windows unicode filepaths are encoded in utf16
        std::string abspath(join_paths(get_cwd(), path));
        if (!(abspath[0] == '/' && abspath[1] == '/')) {
          abspath = "//?/" + abspath;
        }
        std::wstring wpath(UTF_8::convert_to_utf16(abspath));
        std::replace(wpath.begin(), wpath.end(), '/', '\\');
        DWORD rv = GetFullPathNameW(wpath.c_str(), 32767, resolved, NULL);
        if (rv > 32767) throw Exception::OperationError("Path is too long");
        if (rv == 0) throw Exception::OperationError("Path could not be resolved");
        HANDLE hFile = CreateFileW(resolved, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, 0, NULL);
        if (hFile == INVALID_HANDLE_VALUE) return 0;
        DWORD dwFileLength = GetFileSize(hFile, NULL);
        if (dwFileLength == INVALID_FILE_SIZE) return 0;
        // allocate an extra byte for the null char
        // and another one for edge-cases in lexer
        pBuffer = (BYTE*)malloc((dwFileLength+2)*sizeof(BYTE));
        ReadFile(hFile, pBuffer, dwFileLength, &dwBytes, NULL);
        pBuffer[dwFileLength+0] = '\0';
        pBuffer[dwFileLength+1] = '\0';
        CloseHandle(hFile);
        // just convert from unsigned char*
        char* contents = (char*) pBuffer;
      #else
        struct stat st;
        if (stat(path.c_str(), &st) == -1 || S_ISDIR(st.st_mode)) return 0;
        std::ifstream file(path.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
        char* contents = 0;
        if (file.is_open()) {
          size_t size = file.tellg();
          // allocate an extra byte for the null char
          // and another one for edge-cases in lexer
          contents = (char*) malloc((size+2)*sizeof(char));
          file.seekg(0, std::ios::beg);
          file.read(contents, size);
          contents[size+0] = '\0';
          contents[size+1] = '\0';
          file.close();
        }
      #endif
      std::string extension;
      if (path.length() > 5) {
        extension = path.substr(path.length() - 5, 5);
      }
      for(size_t i=0; i<extension.size();++i)
        extension[i] = tolower(extension[i]);
      if (extension == ".sass" && contents != 0) {
        char * converted = sass2scss(contents, SASS2SCSS_PRETTIFY_1 | SASS2SCSS_KEEP_COMMENT);
        free(contents); // free the indented contents
        return converted; // should be freed by caller
      } else {
        return contents;
      }
    }

    // split a path string delimited by semicolons or colons (OS dependent)
    std::vector<std::string> split_path_list(const char* str)
    {
      std::vector<std::string> paths;
      if (str == NULL) return paths;
      // find delimiter via prelexer (return zero at end)
      const char* end = Prelexer::find_first<PATH_SEP>(str);
      // search until null delimiter
      while (end) {
        // add path from current position to delimiter
        paths.push_back(std::string(str, end - str));
        str = end + 1; // skip delimiter
        end = Prelexer::find_first<PATH_SEP>(str);
      }
      // add path from current position to end
      paths.push_back(std::string(str));
      // return back
      return paths;
    }

  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  // simple endless recursion protection
  const size_t maxRecursion = 500;

  Expand::Expand(Context& ctx, Env* env, SelectorStack* stack)
  : ctx(ctx),
    traces(ctx.traces),
    eval(Eval(*this)),
    recursions(0),
    in_keyframes(false),
    at_root_without_rule(false),
    old_at_root_without_rule(false),
    env_stack(EnvStack()),
    block_stack(BlockStack()),
    call_stack(CallStack()),
    selector_stack(SelectorStack()),
    media_stack(MediaStack())
  {
    env_stack.push_back(nullptr);
    env_stack.push_back(env);
    block_stack.push_back(nullptr);
    call_stack.push_back({});
    if (stack == NULL) { selector_stack.push_back({}); }
    else { selector_stack.insert(selector_stack.end(), stack->begin(), stack->end()); }
    media_stack.push_back(nullptr);
  }

  Env* Expand::environment()
  {
    if (env_stack.size() > 0)
      return env_stack.back();
    return 0;
  }

  Selector_List_Obj Expand::selector()
  {
    if (selector_stack.size() > 0)
      return selector_stack.back();
    return {};
  }

  // blocks create new variable scopes
  Block_Ptr Expand::operator()(Block_Ptr b)
  {
    // create new local environment
    // set the current env as parent
    Env env(environment());
    // copy the block object (add items later)
    Block_Obj bb = SASS_MEMORY_NEW(Block,
                                b->pstate(),
                                b->length(),
                                b->is_root());
    // setup block and env stack
    this->block_stack.push_back(bb);
    this->env_stack.push_back(&env);
    // operate on block
    // this may throw up!
    this->append_block(b);
    // revert block and env stack
    this->block_stack.pop_back();
    this->env_stack.pop_back();
    // return copy
    return bb.detach();
  }

  Statement_Ptr Expand::operator()(Ruleset_Ptr r)
  {
    LOCAL_FLAG(old_at_root_without_rule, at_root_without_rule);

    if (in_keyframes) {
      Block_Ptr bb = operator()(r->block());
      Keyframe_Rule_Obj k = SASS_MEMORY_NEW(Keyframe_Rule, r->pstate(), bb);
      if (r->selector()) {
        if (Selector_List_Ptr s = r->selector()) {
          selector_stack.push_back({});
          k->name(s->eval(eval));
          selector_stack.pop_back();
        }
      }
      return k.detach();
    }

    // reset when leaving scope
    LOCAL_FLAG(at_root_without_rule, false);

    // `&` is allowed in `@at-root`!
    bool has_parent_selector = false;
    for (size_t i = 0, L = selector_stack.size(); i < L && !has_parent_selector; i++) {
      Selector_List_Obj ll = selector_stack.at(i);
      has_parent_selector = ll != 0 && ll->length() > 0;
    }

    Selector_List_Obj sel = r->selector();
    if (sel) sel = sel->eval(eval);

    // check for parent selectors in base level rules
    if (r->is_root() || (block_stack.back() && block_stack.back()->is_root())) {
      if (Selector_List_Ptr selector_list = Cast<Selector_List>(r->selector())) {
        for (Complex_Selector_Obj complex_selector : selector_list->elements()) {
          Complex_Selector_Ptr tail = complex_selector;
          while (tail) {
            if (tail->head()) for (Simple_Selector_Obj header : tail->head()->elements()) {
              Parent_Selector_Ptr ptr = Cast<Parent_Selector>(header);
              if (ptr == NULL || (!ptr->real() || has_parent_selector)) continue;
              std::string sel_str(complex_selector->to_string(ctx.c_options));
              error("Base-level rules cannot contain the parent-selector-referencing character '&'.", header->pstate(), traces);
            }
            tail = tail->tail();
          }
        }
      }
    }
    else {
      if (sel->length() == 0 || sel->has_parent_ref()) {
        if (sel->has_real_parent_ref() && !has_parent_selector) {
          error("Base-level rules cannot contain the parent-selector-referencing character '&'.", sel->pstate(), traces);
        }
      }
    }

    // do not connect parent again
    sel->remove_parent_selectors();
    selector_stack.push_back(sel);
    Env env(environment());
    if (block_stack.back()->is_root()) {
      env_stack.push_back(&env);
    }
    sel->set_media_block(media_stack.back());
    Block_Obj blk;
    if (r->block()) blk = operator()(r->block());
    Ruleset_Ptr rr = SASS_MEMORY_NEW(Ruleset,
                                  r->pstate(),
                                  sel,
                                  blk);
    selector_stack.pop_back();
    if (block_stack.back()->is_root()) {
      env_stack.pop_back();
    }

    rr->is_root(r->is_root());
    rr->tabs(r->tabs());

    return rr;
  }

  Statement_Ptr Expand::operator()(Supports_Block_Ptr f)
  {
    Expression_Obj condition = f->condition()->perform(&eval);
    Supports_Block_Obj ff = SASS_MEMORY_NEW(Supports_Block,
                                       f->pstate(),
                                       Cast<Supports_Condition>(condition),
                                       operator()(f->block()));
    return ff.detach();
  }

  Statement_Ptr Expand::operator()(Media_Block_Ptr m)
  {
    Media_Block_Obj cpy = SASS_MEMORY_COPY(m);
    // Media_Blocks are prone to have circular references
    // Copy could leak memory if it does not get picked up
    // Looks like we are able to reset block reference for copy
    // Good as it will ensure a low memory overhead for this fix
    // So this is a cheap solution with a minimal price
    ctx.ast_gc.push_back(cpy); cpy->block({});
    Expression_Obj mq = eval(m->media_queries());
    std::string str_mq(mq->to_string(ctx.c_options));
    char* str = sass_copy_c_string(str_mq.c_str());
    ctx.strings.push_back(str);
    Parser p(Parser::from_c_str(str, ctx, traces, mq->pstate()));
    mq = p.parse_media_queries(); // re-assign now
    cpy->media_queries(mq);
    media_stack.push_back(cpy);
    Block_Obj blk = operator()(m->block());
    Media_Block_Ptr mm = SASS_MEMORY_NEW(Media_Block,
                                      m->pstate(),
                                      mq,
                                      blk);
    media_stack.pop_back();
    mm->tabs(m->tabs());
    return mm;
  }

  Statement_Ptr Expand::operator()(At_Root_Block_Ptr a)
  {
    Block_Obj ab = a->block();
    Expression_Obj ae = a->expression();

    if (ae) ae = ae->perform(&eval);
    else ae = SASS_MEMORY_NEW(At_Root_Query, a->pstate());

    LOCAL_FLAG(at_root_without_rule, Cast<At_Root_Query>(ae)->exclude("rule"));
    LOCAL_FLAG(in_keyframes, false);

                                       ;

    Block_Obj bb = ab ? operator()(ab) : NULL;
    At_Root_Block_Obj aa = SASS_MEMORY_NEW(At_Root_Block,
                                        a->pstate(),
                                        bb,
                                        Cast<At_Root_Query>(ae));
    return aa.detach();
  }

  Statement_Ptr Expand::operator()(Directive_Ptr a)
  {
    LOCAL_FLAG(in_keyframes, a->is_keyframes());
    Block_Ptr ab = a->block();
    Selector_List_Ptr as = a->selector();
    Expression_Ptr av = a->value();
    selector_stack.push_back({});
    if (av) av = av->perform(&eval);
    if (as) as = eval(as);
    selector_stack.pop_back();
    Block_Ptr bb = ab ? operator()(ab) : NULL;
    Directive_Ptr aa = SASS_MEMORY_NEW(Directive,
                                  a->pstate(),
                                  a->keyword(),
                                  as,
                                  bb,
                                  av);
    return aa;
  }

  Statement_Ptr Expand::operator()(Declaration_Ptr d)
  {
    Block_Obj ab = d->block();
    String_Obj old_p = d->property();
    Expression_Obj prop = old_p->perform(&eval);
    String_Obj new_p = Cast<String>(prop);
    // we might get a color back
    if (!new_p) {
      std::string str(prop->to_string(ctx.c_options));
      new_p = SASS_MEMORY_NEW(String_Constant, old_p->pstate(), str);
    }
    Expression_Obj value = d->value();
    if (value) value = value->perform(&eval);
    Block_Obj bb = ab ? operator()(ab) : NULL;
    if (!bb) {
      if (!value || (value->is_invisible() && !d->is_important())) return 0;
    }
    Declaration_Ptr decl = SASS_MEMORY_NEW(Declaration,
                                        d->pstate(),
                                        new_p,
                                        value,
                                        d->is_important(),
                                        d->is_custom_property(),
                                        bb);
    decl->tabs(d->tabs());
    return decl;
  }

  Statement_Ptr Expand::operator()(Assignment_Ptr a)
  {
    Env* env = environment();
    const std::string& var(a->variable());
    if (a->is_global()) {
      if (a->is_default()) {
        if (env->has_global(var)) {
          Expression_Obj e = Cast<Expression>(env->get_global(var));
          if (!e || e->concrete_type() == Expression::NULL_VAL) {
            env->set_global(var, a->value()->perform(&eval));
          }
        }
        else {
          env->set_global(var, a->value()->perform(&eval));
        }
      }
      else {
        env->set_global(var, a->value()->perform(&eval));
      }
    }
    else if (a->is_default()) {
      if (env->has_lexical(var)) {
        auto cur = env;
        while (cur && cur->is_lexical()) {
          if (cur->has_local(var)) {
            if (AST_Node_Obj node = cur->get_local(var)) {
              Expression_Obj e = Cast<Expression>(node);
              if (!e || e->concrete_type() == Expression::NULL_VAL) {
                cur->set_local(var, a->value()->perform(&eval));
              }
            }
            else {
              throw std::runtime_error("Env not in sync");
            }
            return 0;
          }
          cur = cur->parent();
        }
        throw std::runtime_error("Env not in sync");
      }
      else if (env->has_global(var)) {
        if (AST_Node_Obj node = env->get_global(var)) {
          Expression_Obj e = Cast<Expression>(node);
          if (!e || e->concrete_type() == Expression::NULL_VAL) {
            env->set_global(var, a->value()->perform(&eval));
          }
        }
      }
      else if (env->is_lexical()) {
        env->set_local(var, a->value()->perform(&eval));
      }
      else {
        env->set_local(var, a->value()->perform(&eval));
      }
    }
    else {
      env->set_lexical(var, a->value()->perform(&eval));
    }
    return 0;
  }

  Statement_Ptr Expand::operator()(Import_Ptr imp)
  {
    Import_Obj result = SASS_MEMORY_NEW(Import, imp->pstate());
    if (imp->import_queries() && imp->import_queries()->size()) {
      Expression_Obj ex = imp->import_queries()->perform(&eval);
      result->import_queries(Cast<List>(ex));
    }
    for ( size_t i = 0, S = imp->urls().size(); i < S; ++i) {
      result->urls().push_back(imp->urls()[i]->perform(&eval));
    }
    // all resources have been dropped for Input_Stubs
    // for ( size_t i = 0, S = imp->incs().size(); i < S; ++i) {}
    return result.detach();
  }

  Statement_Ptr Expand::operator()(Import_Stub_Ptr i)
  {
    traces.push_back(Backtrace(i->pstate()));
    // get parent node from call stack
    AST_Node_Obj parent = call_stack.back();
    if (Cast<Block>(parent) == NULL) {
      error("Import directives may not be used within control directives or mixins.", i->pstate(), traces);
    }
    // we don't seem to need that actually afterall
    Sass_Import_Entry import = sass_make_import(
      i->imp_path().c_str(),
      i->abs_path().c_str(),
      0, 0
    );
    ctx.import_stack.push_back(import);

    Block_Obj trace_block = SASS_MEMORY_NEW(Block, i->pstate());
    Trace_Obj trace = SASS_MEMORY_NEW(Trace, i->pstate(), i->imp_path(), trace_block, 'i');
    block_stack.back()->append(trace);
    block_stack.push_back(trace_block);

    const std::string& abs_path(i->resource().abs_path);
    append_block(ctx.sheets.at(abs_path).root);
    sass_delete_import(ctx.import_stack.back());
    ctx.import_stack.pop_back();
    block_stack.pop_back();
    traces.pop_back();
    return 0;
  }

  Statement_Ptr Expand::operator()(Warning_Ptr w)
  {
    // eval handles this too, because warnings may occur in functions
    w->perform(&eval);
    return 0;
  }

  Statement_Ptr Expand::operator()(Error_Ptr e)
  {
    // eval handles this too, because errors may occur in functions
    e->perform(&eval);
    return 0;
  }

  Statement_Ptr Expand::operator()(Debug_Ptr d)
  {
    // eval handles this too, because warnings may occur in functions
    d->perform(&eval);
    return 0;
  }

  Statement_Ptr Expand::operator()(Comment_Ptr c)
  {
    if (ctx.output_style() == COMPRESSED) {
      // comments should not be evaluated in compact
      // https://github.com/sass/libsass/issues/2359
      if (!c->is_important()) return NULL;
    }
    eval.is_in_comment = true;
    Comment_Ptr rv = SASS_MEMORY_NEW(Comment, c->pstate(), Cast<String>(c->text()->perform(&eval)), c->is_important());
    eval.is_in_comment = false;
    // TODO: eval the text, once we're parsing/storing it as a String_Schema
    return rv;
  }

  Statement_Ptr Expand::operator()(If_Ptr i)
  {
    Env env(environment(), true);
    env_stack.push_back(&env);
    call_stack.push_back(i);
    Expression_Obj rv = i->predicate()->perform(&eval);
    if (*rv) {
      append_block(i->block());
    }
    else {
      Block_Ptr alt = i->alternative();
      if (alt) append_block(alt);
    }
    call_stack.pop_back();
    env_stack.pop_back();
    return 0;
  }

  // For does not create a new env scope
  // But iteration vars are reset afterwards
  Statement_Ptr Expand::operator()(For_Ptr f)
  {
    std::string variable(f->variable());
    Expression_Obj low = f->lower_bound()->perform(&eval);
    if (low->concrete_type() != Expression::NUMBER) {
      traces.push_back(Backtrace(low->pstate()));
      throw Exception::TypeMismatch(traces, *low, "integer");
    }
    Expression_Obj high = f->upper_bound()->perform(&eval);
    if (high->concrete_type() != Expression::NUMBER) {
      traces.push_back(Backtrace(high->pstate()));
      throw Exception::TypeMismatch(traces, *high, "integer");
    }
    Number_Obj sass_start = Cast<Number>(low);
    Number_Obj sass_end = Cast<Number>(high);
    // check if units are valid for sequence
    if (sass_start->unit() != sass_end->unit()) {
      std::stringstream msg; msg << "Incompatible units: '"
        << sass_start->unit() << "' and '"
        << sass_end->unit() << "'.";
      error(msg.str(), low->pstate(), traces);
    }
    double start = sass_start->value();
    double end = sass_end->value();
    // only create iterator once in this environment
    Env env(environment(), true);
    env_stack.push_back(&env);
    call_stack.push_back(f);
    Block_Ptr body = f->block();
    if (start < end) {
      if (f->is_inclusive()) ++end;
      for (double i = start;
           i < end;
           ++i) {
        Number_Obj it = SASS_MEMORY_NEW(Number, low->pstate(), i, sass_end->unit());
        env.set_local(variable, it);
        append_block(body);
      }
    } else {
      if (f->is_inclusive()) --end;
      for (double i = start;
           i > end;
           --i) {
        Number_Obj it = SASS_MEMORY_NEW(Number, low->pstate(), i, sass_end->unit());
        env.set_local(variable, it);
        append_block(body);
      }
    }
    call_stack.pop_back();
    env_stack.pop_back();
    return 0;
  }

  // Eval does not create a new env scope
  // But iteration vars are reset afterwards
  Statement_Ptr Expand::operator()(Each_Ptr e)
  {
    std::vector<std::string> variables(e->variables());
    Expression_Obj expr = e->list()->perform(&eval);
    List_Obj list;
    Map_Obj map;
    if (expr->concrete_type() == Expression::MAP) {
      map = Cast<Map>(expr);
    }
    else if (Selector_List_Ptr ls = Cast<Selector_List>(expr)) {
      Listize listize;
      Expression_Obj rv = ls->perform(&listize);
      list = Cast<List>(rv);
    }
    else if (expr->concrete_type() != Expression::LIST) {
      list = SASS_MEMORY_NEW(List, expr->pstate(), 1, SASS_COMMA);
      list->append(expr);
    }
    else {
      list = Cast<List>(expr);
    }
    // remember variables and then reset them
    Env env(environment(), true);
    env_stack.push_back(&env);
    call_stack.push_back(e);
    Block_Ptr body = e->block();

    if (map) {
      for (auto key : map->keys()) {
        Expression_Obj k = key->perform(&eval);
        Expression_Obj v = map->at(key)->perform(&eval);

        if (variables.size() == 1) {
          List_Obj variable = SASS_MEMORY_NEW(List, map->pstate(), 2, SASS_SPACE);
          variable->append(k);
          variable->append(v);
          env.set_local(variables[0], variable);
        } else {
          env.set_local(variables[0], k);
          env.set_local(variables[1], v);
        }
        append_block(body);
      }
    }
    else {
      // bool arglist = list->is_arglist();
      if (list->length() == 1 && Cast<Selector_List>(list)) {
        list = Cast<List>(list);
      }
      for (size_t i = 0, L = list->length(); i < L; ++i) {
        Expression_Obj item = list->at(i);
        // unwrap value if the expression is an argument
        if (Argument_Obj arg = Cast<Argument>(item)) item = arg->value();
        // check if we got passed a list of args (investigate)
        if (List_Obj scalars = Cast<List>(item)) {
          if (variables.size() == 1) {
            List_Obj var = scalars;
            // if (arglist) var = (*scalars)[0];
            env.set_local(variables[0], var);
          } else {
            for (size_t j = 0, K = variables.size(); j < K; ++j) {
              Expression_Obj res = j >= scalars->length()
                ? SASS_MEMORY_NEW(Null, expr->pstate())
                : (*scalars)[j]->perform(&eval);
              env.set_local(variables[j], res);
            }
          }
        } else {
          if (variables.size() > 0) {
            env.set_local(variables.at(0), item);
            for (size_t j = 1, K = variables.size(); j < K; ++j) {
              Expression_Obj res = SASS_MEMORY_NEW(Null, expr->pstate());
              env.set_local(variables[j], res);
            }
          }
        }
        append_block(body);
      }
    }
    call_stack.pop_back();
    env_stack.pop_back();
    return 0;
  }

  Statement_Ptr Expand::operator()(While_Ptr w)
  {
    Expression_Obj pred = w->predicate();
    Block_Ptr body = w->block();
    Env env(environment(), true);
    env_stack.push_back(&env);
    call_stack.push_back(w);
    Expression_Obj cond = pred->perform(&eval);
    while (!cond->is_false()) {
      append_block(body);
      cond = pred->perform(&eval);
    }
    call_stack.pop_back();
    env_stack.pop_back();
    return 0;
  }

  Statement_Ptr Expand::operator()(Return_Ptr r)
  {
    error("@return may only be used within a function", r->pstate(), traces);
    return 0;
  }


  void Expand::expand_selector_list(Selector_Obj s, Selector_List_Obj extender) {

    if (Selector_List_Obj sl = Cast<Selector_List>(s)) {
      for (Complex_Selector_Obj complex_selector : sl->elements()) {
        Complex_Selector_Obj tail = complex_selector;
        while (tail) {
          if (tail->head()) for (Simple_Selector_Obj header : tail->head()->elements()) {
            if (Cast<Parent_Selector>(header) == NULL) continue; // skip all others
            std::string sel_str(complex_selector->to_string(ctx.c_options));
            error("Can't extend " + sel_str + ": can't extend parent selectors", header->pstate(), traces);
          }
          tail = tail->tail();
        }
      }
    }


    Selector_List_Obj contextualized = Cast<Selector_List>(s->perform(&eval));
    if (contextualized == false) return;
    for (auto complex_sel : contextualized->elements()) {
      Complex_Selector_Obj c = complex_sel;
      if (!c->head() || c->tail()) {
        std::string sel_str(contextualized->to_string(ctx.c_options));
        error("Can't extend " + sel_str + ": can't extend nested selectors", c->pstate(), traces);
      }
      Compound_Selector_Obj target = c->head();
      if (contextualized->is_optional()) target->is_optional(true);
      for (size_t i = 0, L = extender->length(); i < L; ++i) {
        Complex_Selector_Obj sel = (*extender)[i];
        if (!(sel->head() && sel->head()->length() > 0 &&
            Cast<Parent_Selector>((*sel->head())[0])))
        {
          Compound_Selector_Obj hh = SASS_MEMORY_NEW(Compound_Selector, (*extender)[i]->pstate());
          hh->media_block((*extender)[i]->media_block());
          Complex_Selector_Obj ssel = SASS_MEMORY_NEW(Complex_Selector, (*extender)[i]->pstate());
          ssel->media_block((*extender)[i]->media_block());
          if (sel->has_line_feed()) ssel->has_line_feed(true);
          Parent_Selector_Obj ps = SASS_MEMORY_NEW(Parent_Selector, (*extender)[i]->pstate());
          ps->media_block((*extender)[i]->media_block());
          hh->append(ps);
          ssel->tail(sel);
          ssel->head(hh);
          sel = ssel;
        }
        // if (c->has_line_feed()) sel->has_line_feed(true);
        ctx.subset_map.put(target, std::make_pair(sel, target));
      }
    }

  }

  Statement* Expand::operator()(Extension_Ptr e)
  {
    if (Selector_List_Obj extender = selector()) {
      Selector_List_Ptr sl = e->selector();
      // abort on invalid selector
      if (sl == NULL) return NULL;
      if (Selector_Schema_Ptr schema = sl->schema()) {
        if (schema->has_real_parent_ref()) {
          // put root block on stack again (ignore parents)
          // selector schema must not connect in eval!
          block_stack.push_back(block_stack.at(1));
          sl = eval(sl->schema());
          block_stack.pop_back();
        } else {
          selector_stack.push_back({});
          sl = eval(sl->schema());
          selector_stack.pop_back();
        }
      }
      for (Complex_Selector_Obj cs : sl->elements()) {
        if (!cs.isNull() && !cs->head().isNull()) {
          cs->head()->media_block(media_stack.back());
        }
      }
      selector_stack.push_back({});
      expand_selector_list(sl, extender);
      selector_stack.pop_back();
    }
    return 0;
  }

  Statement_Ptr Expand::operator()(Definition_Ptr d)
  {
    Env* env = environment();
    Definition_Obj dd = SASS_MEMORY_COPY(d);
    env->local_frame()[d->name() +
                        (d->type() == Definition::MIXIN ? "[m]" : "[f]")] = dd;

    if (d->type() == Definition::FUNCTION && (
      Prelexer::calc_fn_call(d->name().c_str()) ||
      d->name() == "element"    ||
      d->name() == "expression" ||
      d->name() == "url"
    )) {
      deprecated(
        "Naming a function \"" + d->name() + "\" is disallowed and will be an error in future versions of Sass.",
        "This name conflicts with an existing CSS function with special parse rules.",
        false, d->pstate()
      );
    }

    // set the static link so we can have lexical scoping
    dd->environment(env);
    return 0;
  }

  Statement_Ptr Expand::operator()(Mixin_Call_Ptr c)
  {
    if (recursions > maxRecursion) {
      throw Exception::StackError(traces, *c);
    }

    recursions ++;

    Env* env = environment();
    std::string full_name(c->name() + "[m]");
    if (!env->has(full_name)) {
      error("no mixin named " + c->name(), c->pstate(), traces);
    }
    Definition_Obj def = Cast<Definition>((*env)[full_name]);
    Block_Obj body = def->block();
    Parameters_Obj params = def->parameters();

    if (c->block() && c->name() != "@content" && !body->has_content()) {
      error("Mixin \"" + c->name() + "\" does not accept a content block.", c->pstate(), traces);
    }
    Expression_Obj rv = c->arguments()->perform(&eval);
    Arguments_Obj args = Cast<Arguments>(rv);
    std::string msg(", in mixin `" + c->name() + "`");
    traces.push_back(Backtrace(c->pstate(), msg));
    ctx.callee_stack.push_back({
      c->name().c_str(),
      c->pstate().path,
      c->pstate().line + 1,
      c->pstate().column + 1,
      SASS_CALLEE_MIXIN,
      { env }
    });

    Env new_env(def->environment());
    env_stack.push_back(&new_env);
    if (c->block()) {
      Parameters_Obj params = c->block_parameters();
      if (!params) params = SASS_MEMORY_NEW(Parameters, c->pstate());
      // represent mixin content blocks as thunks/closures
      Definition_Obj thunk = SASS_MEMORY_NEW(Definition,
                                          c->pstate(),
                                          "@content",
                                          params,
                                          c->block(),
                                          Definition::MIXIN);
      thunk->environment(env);
      new_env.local_frame()["@content[m]"] = thunk;
    }

    bind(std::string("Mixin"), c->name(), params, args, &new_env, &eval, traces);

    Block_Obj trace_block = SASS_MEMORY_NEW(Block, c->pstate());
    Trace_Obj trace = SASS_MEMORY_NEW(Trace, c->pstate(), c->name(), trace_block);

    env->set_global("is_in_mixin", bool_true);
    if (Block_Ptr pr = block_stack.back()) {
      trace_block->is_root(pr->is_root());
    }
    block_stack.push_back(trace_block);
    for (auto bb : body->elements()) {
      if (Ruleset_Ptr r = Cast<Ruleset>(bb)) {
        r->is_root(trace_block->is_root());
      }
      Statement_Obj ith = bb->perform(this);
      if (ith) trace->block()->append(ith);
    }
    block_stack.pop_back();
    env->del_global("is_in_mixin");

    ctx.callee_stack.pop_back();
    env_stack.pop_back();
    traces.pop_back();

    recursions --;
    return trace.detach();
  }

  Statement_Ptr Expand::operator()(Content_Ptr c)
  {
    Env* env = environment();
    // convert @content directives into mixin calls to the underlying thunk
    if (!env->has("@content[m]")) return 0;

    if (block_stack.back()->is_root()) {
      selector_stack.push_back({});
    }

    Arguments_Obj args = c->arguments();
    if (!args) args = SASS_MEMORY_NEW(Arguments, c->pstate());

    Mixin_Call_Obj call = SASS_MEMORY_NEW(Mixin_Call,
                                       c->pstate(),
                                       "@content",
                                       args);

    Trace_Obj trace = Cast<Trace>(call->perform(this));

    if (block_stack.back()->is_root()) {
      selector_stack.pop_back();
    }

    return trace.detach();
  }

  // process and add to last block on stack
  inline void Expand::append_block(Block_Ptr b)
  {
    if (b->is_root()) call_stack.push_back(b);
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Ptr stm = b->at(i);
      Statement_Obj ith = stm->perform(this);
      if (ith) block_stack.back()->append(ith);
    }
    if (b->is_root()) call_stack.pop_back();
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {
  namespace UTF_8 {
    using std::string;

    // naming conventions:
    // offset: raw byte offset (0 based)
    // position: code point offset (0 based)
    // index: code point offset (1 based or negative)

    // function that will count the number of code points (utf-8 characters) from the given beginning to the given end
    size_t code_point_count(const string& str, size_t start, size_t end) {
      return utf8::distance(str.begin() + start, str.begin() + end);
    }

    size_t code_point_count(const string& str) {
      return utf8::distance(str.begin(), str.end());
    }

    // function that will return the byte offset at a code point position
    size_t offset_at_position(const string& str, size_t position) {
      string::const_iterator it = str.begin();
      utf8::advance(it, position, str.end());
      return std::distance(str.begin(), it);
    }

    // function that returns number of bytes in a character at offset
    size_t code_point_size_at_offset(const string& str, size_t offset) {
      // get iterator from string and forward by offset
      string::const_iterator stop = str.begin() + offset;
      // check if beyond boundary
      if (stop == str.end()) return 0;
      // advance by one code point
      utf8::advance(stop, 1, str.end());
      // calculate offset for code point
      return  stop - str.begin() - offset;
    }

    // function that will return a normalized index, given a crazy one
    size_t normalize_index(int index, size_t len) {
      long signed_len = static_cast<long>(len);
      // assuming the index is 1-based
      // we are returning a 0-based index
      if (index > 0 && index <= signed_len) {
        // positive and within string length
        return index-1;
      }
      else if (index > signed_len) {
        // positive and past string length
        return len;
      }
      else if (index == 0) {
        return 0;
      }
      else if (std::abs((double)index) <= signed_len) {
        // negative and within string length
        return index + signed_len;
      }
      else {
        // negative and past string length
        return 0;
      }
    }

    #ifdef _WIN32

    // utf16 functions
    using std::wstring;

    // convert from utf16/wide string to utf8 string
    string convert_from_utf16(const wstring& utf16)
    {
      string utf8;
      // pre-allocate expected memory
      utf8.reserve(sizeof(utf16)/2);
      utf8::utf16to8(utf16.begin(), utf16.end(),
                     back_inserter(utf8));
      return utf8;
    }

    // convert from utf8 string to utf16/wide string
    wstring convert_to_utf16(const string& utf8)
    {
      wstring utf16;
      // pre-allocate expected memory
      utf16.reserve(code_point_count(utf8)*2);
      utf8::utf8to16(utf8.begin(), utf8.end(),
                     back_inserter(utf16));
      return utf16;
    }

    #endif

  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  namespace Functions {

    void handle_utf8_error (const ParserState& pstate, Backtraces traces)
    {
      try {
       throw;
      }
      catch (utf8::invalid_code_point&) {
        std::string msg("utf8::invalid_code_point");
        error(msg, pstate, traces);
      }
      catch (utf8::not_enough_room&) {
        std::string msg("utf8::not_enough_room");
        error(msg, pstate, traces);
      }
      catch (utf8::invalid_utf8&) {
        std::string msg("utf8::invalid_utf8");
        error(msg, pstate, traces);
      }
      catch (...) { throw; }
    }

    ///////////////////
    // STRING FUNCTIONS
    ///////////////////

    Signature unquote_sig = "unquote($string)";
    BUILT_IN(sass_unquote)
    {
      AST_Node_Obj arg = env["$string"];
      if (String_Quoted_Ptr string_quoted = Cast<String_Quoted>(arg)) {
        String_Constant_Ptr result = SASS_MEMORY_NEW(String_Constant, pstate, string_quoted->value());
        // remember if the string was quoted (color tokens)
        result->is_delayed(true); // delay colors
        return result;
      }
      else if (String_Constant_Ptr str = Cast<String_Constant>(arg)) {
        return str;
      }
      else if (Value_Ptr ex = Cast<Value>(arg)) {
        Sass_Output_Style oldstyle = ctx.c_options.output_style;
        ctx.c_options.output_style = SASS_STYLE_NESTED;
        std::string val(arg->to_string(ctx.c_options));
        val = Cast<Null>(arg) ? "null" : val;
        ctx.c_options.output_style = oldstyle;

        deprecated_function("Passing " + val + ", a non-string value, to unquote()", pstate);
        return ex;
      }
      throw std::runtime_error("Invalid Data Type for unquote");
    }

    Signature quote_sig = "quote($string)";
    BUILT_IN(sass_quote)
    {
      AST_Node_Obj arg = env["$string"];
      // only set quote mark to true if already a string
      if (String_Quoted_Ptr qstr = Cast<String_Quoted>(arg)) {
        qstr->quote_mark('*');
        return qstr;
      }
      // all other nodes must be converted to a string node
      std::string str(quote(arg->to_string(ctx.c_options), '"'));
      String_Quoted_Ptr result = SASS_MEMORY_NEW(String_Quoted, pstate, str);
      result->quote_mark('*');
      return result;
    }


    Signature str_length_sig = "str-length($string)";
    BUILT_IN(str_length)
    {
      size_t len = std::string::npos;
      try {
        String_Constant_Ptr s = ARG("$string", String_Constant);
        len = UTF_8::code_point_count(s->value(), 0, s->value().size());

      }
      // handle any invalid utf8 errors
      // other errors will be re-thrown
      catch (...) { handle_utf8_error(pstate, traces); }
      // return something even if we had an error (-1)
      return SASS_MEMORY_NEW(Number, pstate, (double)len);
    }

    Signature str_insert_sig = "str-insert($string, $insert, $index)";
    BUILT_IN(str_insert)
    {
      std::string str;
      try {
        String_Constant_Ptr s = ARG("$string", String_Constant);
        str = s->value();
        str = unquote(str);
        String_Constant_Ptr i = ARG("$insert", String_Constant);
        std::string ins = i->value();
        ins = unquote(ins);
        double index = ARGVAL("$index");
        size_t len = UTF_8::code_point_count(str, 0, str.size());

        if (index > 0 && index <= len) {
          // positive and within string length
          str.insert(UTF_8::offset_at_position(str, static_cast<size_t>(index) - 1), ins);
        }
        else if (index > len) {
          // positive and past string length
          str += ins;
        }
        else if (index == 0) {
          str = ins + str;
        }
        else if (std::abs(index) <= len) {
          // negative and within string length
          index += len + 1;
          str.insert(UTF_8::offset_at_position(str, static_cast<size_t>(index)), ins);
        }
        else {
          // negative and past string length
          str = ins + str;
        }

        if (String_Quoted_Ptr ss = Cast<String_Quoted>(s)) {
          if (ss->quote_mark()) str = quote(str);
        }
      }
      // handle any invalid utf8 errors
      // other errors will be re-thrown
      catch (...) { handle_utf8_error(pstate, traces); }
      return SASS_MEMORY_NEW(String_Quoted, pstate, str);
    }

    Signature str_index_sig = "str-index($string, $substring)";
    BUILT_IN(str_index)
    {
      size_t index = std::string::npos;
      try {
        String_Constant_Ptr s = ARG("$string", String_Constant);
        String_Constant_Ptr t = ARG("$substring", String_Constant);
        std::string str = s->value();
        str = unquote(str);
        std::string substr = t->value();
        substr = unquote(substr);

        size_t c_index = str.find(substr);
        if(c_index == std::string::npos) {
          return SASS_MEMORY_NEW(Null, pstate);
        }
        index = UTF_8::code_point_count(str, 0, c_index) + 1;
      }
      // handle any invalid utf8 errors
      // other errors will be re-thrown
      catch (...) { handle_utf8_error(pstate, traces); }
      // return something even if we had an error (-1)
      return SASS_MEMORY_NEW(Number, pstate, (double)index);
    }

    Signature str_slice_sig = "str-slice($string, $start-at, $end-at:-1)";
    BUILT_IN(str_slice)
    {
      std::string newstr;
      try {
        String_Constant_Ptr s = ARG("$string", String_Constant);
        double start_at = ARGVAL("$start-at");
        double end_at = ARGVAL("$end-at");
        String_Quoted_Ptr ss = Cast<String_Quoted>(s);

        std::string str = unquote(s->value());

        size_t size = utf8::distance(str.begin(), str.end());

        if (!Cast<Number>(env["$end-at"])) {
          end_at = -1;
        }

        if (end_at == 0 || (end_at + size) < 0) {
          if (ss && ss->quote_mark()) newstr = quote("");
          return SASS_MEMORY_NEW(String_Quoted, pstate, newstr);
        }

        if (end_at < 0) {
          end_at += size + 1;
          if (end_at == 0) end_at = 1;
        }
        if (end_at > size) { end_at = (double)size; }
        if (start_at < 0) {
          start_at += size + 1;
          if (start_at < 0)  start_at = 0;
        }
        else if (start_at == 0) { ++ start_at; }

        if (start_at <= end_at)
        {
          std::string::iterator start = str.begin();
          utf8::advance(start, start_at - 1, str.end());
          std::string::iterator end = start;
          utf8::advance(end, end_at - start_at + 1, str.end());
          newstr = std::string(start, end);
        }
        if (ss) {
          if(ss->quote_mark()) newstr = quote(newstr);
        }
      }
      // handle any invalid utf8 errors
      // other errors will be re-thrown
      catch (...) { handle_utf8_error(pstate, traces); }
      return SASS_MEMORY_NEW(String_Quoted, pstate, newstr);
    }

    Signature to_upper_case_sig = "to-upper-case($string)";
    BUILT_IN(to_upper_case)
    {
      String_Constant_Ptr s = ARG("$string", String_Constant);
      std::string str = s->value();

      for (size_t i = 0, L = str.length(); i < L; ++i) {
        if (Sass::Util::isAscii(str[i])) {
          str[i] = std::toupper(str[i]);
        }
      }

      if (String_Quoted_Ptr ss = Cast<String_Quoted>(s)) {
        String_Quoted_Ptr cpy = SASS_MEMORY_COPY(ss);
        cpy->value(str);
        return cpy;
      } else {
        return SASS_MEMORY_NEW(String_Quoted, pstate, str);
      }
    }

    Signature to_lower_case_sig = "to-lower-case($string)";
    BUILT_IN(to_lower_case)
    {
      String_Constant_Ptr s = ARG("$string", String_Constant);
      std::string str = s->value();

      for (size_t i = 0, L = str.length(); i < L; ++i) {
        if (Sass::Util::isAscii(str[i])) {
          str[i] = std::tolower(str[i]);
        }
      }

      if (String_Quoted_Ptr ss = Cast<String_Quoted>(s)) {
        String_Quoted_Ptr cpy = SASS_MEMORY_COPY(ss);
        cpy->value(str);
        return cpy;
      } else {
        return SASS_MEMORY_NEW(String_Quoted, pstate, str);
      }
    }

  }

}
/*
  Copyright (C) 2011 Joseph A. Adams (joeyadams3.14159@gmail.com)
  All rights reserved.

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.
*/

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_DEPRECATE
#endif


// include utf8 library used by libsass
// ToDo: replace internal json utf8 code

#include <string.h>

#if defined(_MSC_VER) && _MSC_VER < 1900
#ifdef snprintf
#undef snprintf
#endif
extern "C" int snprintf(char *, size_t, const char *, ...);
#endif

#define out_of_memory() do {                    \
    fprintf(stderr, "Out of memory.\n");    \
    exit(EXIT_FAILURE);                     \
  } while (0)

/* Sadly, strdup is not portable. */
static char *json_strdup(const char *str)
{
  char *ret = (char*) malloc(strlen(str) + 1);
  if (ret == NULL)
    out_of_memory();
  strcpy(ret, str);
  return ret;
}

/* String buffer */

typedef struct
{
  char *cur;
  char *end;
  char *start;
} SB;

static void sb_init(SB *sb)
{
  sb->start = (char*) malloc(17);
  if (sb->start == NULL)
    out_of_memory();
  sb->cur = sb->start;
  sb->end = sb->start + 16;
}

/* sb and need may be evaluated multiple times. */
#define sb_need(sb, need) do {                  \
    if ((sb)->end - (sb)->cur < (need))     \
      sb_grow(sb, need);                  \
  } while (0)

static void sb_grow(SB *sb, int need)
{
  size_t length = sb->cur - sb->start;
  size_t alloc = sb->end - sb->start;

  do {
    alloc *= 2;
  } while (alloc < length + need);

  sb->start = (char*) realloc(sb->start, alloc + 1);
  if (sb->start == NULL)
    out_of_memory();
  sb->cur = sb->start + length;
  sb->end = sb->start + alloc;
}

static void sb_put(SB *sb, const char *bytes, int count)
{
  sb_need(sb, count);
  memcpy(sb->cur, bytes, count);
  sb->cur += count;
}

#define sb_putc(sb, c) do {         \
    if ((sb)->cur >= (sb)->end) \
      sb_grow(sb, 1);         \
    *(sb)->cur++ = (c);         \
  } while (0)

static void sb_puts(SB *sb, const char *str)
{
  sb_put(sb, str, (int)strlen(str));
}

static char *sb_finish(SB *sb)
{
  *sb->cur = 0;
  assert(sb->start <= sb->cur && strlen(sb->start) == (size_t)(sb->cur - sb->start));
  return sb->start;
}

static void sb_free(SB *sb)
{
  free(sb->start);
}

/*
 * Unicode helper functions
 *
 * These are taken from the ccan/charset module and customized a bit.
 * Putting them here means the compiler can (choose to) inline them,
 * and it keeps ccan/json from having a dependency.
 *
 * We use uint32_t Type for Unicode codepoints.
 * We need our own because wchar_t might be 16 bits.
 */

/*
 * Validate a single UTF-8 character starting at @s.
 * The string must be null-terminated.
 *
 * If it's valid, return its length (1 thru 4).
 * If it's invalid or clipped, return 0.
 *
 * This function implements the syntax given in RFC3629, which is
 * the same as that given in The Unicode Standard, Version 6.0.
 *
 * It has the following properties:
 *
 *  * All codepoints U+0000..U+10FFFF may be encoded,
 *    except for U+D800..U+DFFF, which are reserved
 *    for UTF-16 surrogate pair encoding.
 *  * UTF-8 byte sequences longer than 4 bytes are not permitted,
 *    as they exceed the range of Unicode.
 *  * The sixty-six Unicode "non-characters" are permitted
 *    (namely, U+FDD0..U+FDEF, U+xxFFFE, and U+xxFFFF).
 */
static int utf8_validate_cz(const char *s)
{
  unsigned char c = *s++;

  if (c <= 0x7F) {        /* 00..7F */
    return 1;
  } else if (c <= 0xC1) { /* 80..C1 */
    /* Disallow overlong 2-byte sequence. */
    return 0;
  } else if (c <= 0xDF) { /* C2..DF */
    /* Make sure subsequent byte is in the range 0x80..0xBF. */
    if (((unsigned char)*s++ & 0xC0) != 0x80)
      return 0;

    return 2;
  } else if (c <= 0xEF) { /* E0..EF */
    /* Disallow overlong 3-byte sequence. */
    if (c == 0xE0 && (unsigned char)*s < 0xA0)
      return 0;

    /* Disallow U+D800..U+DFFF. */
    if (c == 0xED && (unsigned char)*s > 0x9F)
      return 0;

    /* Make sure subsequent bytes are in the range 0x80..0xBF. */
    if (((unsigned char)*s++ & 0xC0) != 0x80)
      return 0;
    if (((unsigned char)*s++ & 0xC0) != 0x80)
      return 0;

    return 3;
  } else if (c <= 0xF4) { /* F0..F4 */
    /* Disallow overlong 4-byte sequence. */
    if (c == 0xF0 && (unsigned char)*s < 0x90)
      return 0;

    /* Disallow codepoints beyond U+10FFFF. */
    if (c == 0xF4 && (unsigned char)*s > 0x8F)
      return 0;

    /* Make sure subsequent bytes are in the range 0x80..0xBF. */
    if (((unsigned char)*s++ & 0xC0) != 0x80)
      return 0;
    if (((unsigned char)*s++ & 0xC0) != 0x80)
      return 0;
    if (((unsigned char)*s++ & 0xC0) != 0x80)
      return 0;

    return 4;
  } else {                /* F5..FF */
    return 0;
  }
}

/* Validate a null-terminated UTF-8 string. */
static bool utf8_validate(const char *s)
{
  int len;

  for (; *s != 0; s += len) {
    len = utf8_validate_cz(s);
    if (len == 0)
      return false;
  }

  return true;
}

/*
 * Read a single UTF-8 character starting at @s,
 * returning the length, in bytes, of the character read.
 *
 * This function assumes input is valid UTF-8,
 * and that there are enough characters in front of @s.
 */
static int utf8_read_char(const char *s, uint32_t *out)
{
  const unsigned char *c = (const unsigned char*) s;

  assert(utf8_validate_cz(s));

  if (c[0] <= 0x7F) {
    /* 00..7F */
    *out = c[0];
    return 1;
  } else if (c[0] <= 0xDF) {
    /* C2..DF (unless input is invalid) */
    *out = ((uint32_t)c[0] & 0x1F) << 6 |
           ((uint32_t)c[1] & 0x3F);
    return 2;
  } else if (c[0] <= 0xEF) {
    /* E0..EF */
    *out = ((uint32_t)c[0] &  0xF) << 12 |
           ((uint32_t)c[1] & 0x3F) << 6  |
           ((uint32_t)c[2] & 0x3F);
    return 3;
  } else {
    /* F0..F4 (unless input is invalid) */
    *out = ((uint32_t)c[0] &  0x7) << 18 |
           ((uint32_t)c[1] & 0x3F) << 12 |
           ((uint32_t)c[2] & 0x3F) << 6  |
           ((uint32_t)c[3] & 0x3F);
    return 4;
  }
}

/*
 * Write a single UTF-8 character to @s,
 * returning the length, in bytes, of the character written.
 *
 * @unicode must be U+0000..U+10FFFF, but not U+D800..U+DFFF.
 *
 * This function will write up to 4 bytes to @out.
 */
static int utf8_write_char(uint32_t unicode, char *out)
{
  unsigned char *o = (unsigned char*) out;

  assert(unicode <= 0x10FFFF && !(unicode >= 0xD800 && unicode <= 0xDFFF));

  if (unicode <= 0x7F) {
    /* U+0000..U+007F */
    *o++ = unicode;
    return 1;
  } else if (unicode <= 0x7FF) {
    /* U+0080..U+07FF */
    *o++ = 0xC0 | unicode >> 6;
    *o++ = 0x80 | (unicode & 0x3F);
    return 2;
  } else if (unicode <= 0xFFFF) {
    /* U+0800..U+FFFF */
    *o++ = 0xE0 | unicode >> 12;
    *o++ = 0x80 | (unicode >> 6 & 0x3F);
    *o++ = 0x80 | (unicode & 0x3F);
    return 3;
  } else {
    /* U+10000..U+10FFFF */
    *o++ = 0xF0 | unicode >> 18;
    *o++ = 0x80 | (unicode >> 12 & 0x3F);
    *o++ = 0x80 | (unicode >> 6 & 0x3F);
    *o++ = 0x80 | (unicode & 0x3F);
    return 4;
  }
}

/*
 * Compute the Unicode codepoint of a UTF-16 surrogate pair.
 *
 * @uc should be 0xD800..0xDBFF, and @lc should be 0xDC00..0xDFFF.
 * If they aren't, this function returns false.
 */
static bool from_surrogate_pair(uint16_t uc, uint16_t lc, uint32_t *unicode)
{
  if (uc >= 0xD800 && uc <= 0xDBFF && lc >= 0xDC00 && lc <= 0xDFFF) {
    *unicode = 0x10000 + ((((uint32_t)uc & 0x3FF) << 10) | (lc & 0x3FF));
    return true;
  } else {
    return false;
  }
}

/*
 * Construct a UTF-16 surrogate pair given a Unicode codepoint.
 *
 * @unicode must be U+10000..U+10FFFF.
 */
static void to_surrogate_pair(uint32_t unicode, uint16_t *uc, uint16_t *lc)
{
  uint32_t n;

  assert(unicode >= 0x10000 && unicode <= 0x10FFFF);

  n = unicode - 0x10000;
  *uc = ((n >> 10) & 0x3FF) | 0xD800;
  *lc = (n & 0x3FF) | 0xDC00;
}

static bool is_space        (const char *c);
static bool is_digit        (const char *c);
static bool parse_value     (const char **sp, JsonNode        **out);
static bool parse_string    (const char **sp, char            **out);
static bool parse_number    (const char **sp, double           *out);
static bool parse_array     (const char **sp, JsonNode        **out);
static bool parse_object    (const char **sp, JsonNode        **out);
static bool parse_hex16     (const char **sp, uint16_t         *out);

static bool expect_literal  (const char **sp, const char *str);
static void skip_space      (const char **sp);

static void emit_value              (SB *out, const JsonNode *node);
static void emit_value_indented     (SB *out, const JsonNode *node, const char *space, int indent_level);
static void emit_string             (SB *out, const char *str);
static void emit_number             (SB *out, double num);
static void emit_array              (SB *out, const JsonNode *array);
static void emit_array_indented     (SB *out, const JsonNode *array, const char *space, int indent_level);
static void emit_object             (SB *out, const JsonNode *object);
static void emit_object_indented    (SB *out, const JsonNode *object, const char *space, int indent_level);

static int write_hex16(char *out, uint16_t val);

static JsonNode *mknode(JsonTag tag);
static void append_node(JsonNode *parent, JsonNode *child);
static void prepend_node(JsonNode *parent, JsonNode *child);
static void append_member(JsonNode *object, char *key, JsonNode *value);

/* Assertion-friendly validity checks */
static bool tag_is_valid(unsigned int tag);
static bool number_is_valid(const char *num);

JsonNode *json_decode(const char *json)
{
  const char *s = json;
  JsonNode *ret;

  skip_space(&s);
  if (!parse_value(&s, &ret))
    return NULL;

  skip_space(&s);
  if (*s != 0) {
    json_delete(ret);
    return NULL;
  }

  return ret;
}

char *json_encode(const JsonNode *node)
{
  return json_stringify(node, NULL);
}

char *json_encode_string(const char *str)
{
  SB sb;
  sb_init(&sb);

  try {
    emit_string(&sb, str);
  }
  catch (std::exception&) {
    sb_free(&sb);
    throw;
  }

  return sb_finish(&sb);
}

char *json_stringify(const JsonNode *node, const char *space)
{
  SB sb;
  sb_init(&sb);

  try {
    if (space != NULL)
      emit_value_indented(&sb, node, space, 0);
    else
      emit_value(&sb, node);
  }
  catch (std::exception&) {
    sb_free(&sb);
    throw;
  }

  return sb_finish(&sb);
}

void json_delete(JsonNode *node)
{
  if (node != NULL) {
    json_remove_from_parent(node);

    switch (node->tag) {
      case JSON_STRING:
        free(node->string_);
        break;
      case JSON_ARRAY:
      case JSON_OBJECT:
      {
        JsonNode *child, *next;
        for (child = node->children.head; child != NULL; child = next) {
          next = child->next;
          json_delete(child);
        }
        break;
      }
      default:;
    }

    free(node);
  }
}

bool json_validate(const char *json)
{
  const char *s = json;

  skip_space(&s);
  if (!parse_value(&s, NULL))
    return false;

  skip_space(&s);
  if (*s != 0)
    return false;

  return true;
}

JsonNode *json_find_element(JsonNode *array, int index)
{
  JsonNode *element;
  int i = 0;

  if (array == NULL || array->tag != JSON_ARRAY)
    return NULL;

  json_foreach(element, array) {
    if (i == index)
      return element;
    i++;
  }

  return NULL;
}

JsonNode *json_find_member(JsonNode *object, const char *name)
{
  JsonNode *member;

  if (object == NULL || object->tag != JSON_OBJECT)
    return NULL;

  json_foreach(member, object)
    if (strcmp(member->key, name) == 0)
      return member;

  return NULL;
}

JsonNode *json_first_child(const JsonNode *node)
{
  if (node != NULL && (node->tag == JSON_ARRAY || node->tag == JSON_OBJECT))
    return node->children.head;
  return NULL;
}

static JsonNode *mknode(JsonTag tag)
{
  JsonNode *ret = (JsonNode*) calloc(1, sizeof(JsonNode));
  if (ret == NULL)
    out_of_memory();
  ret->tag = tag;
  return ret;
}

JsonNode *json_mknull(void)
{
  return mknode(JSON_NULL);
}

JsonNode *json_mkbool(bool b)
{
  JsonNode *ret = mknode(JSON_BOOL);
  ret->bool_ = b;
  return ret;
}

static JsonNode *mkstring(char *s)
{
  JsonNode *ret = mknode(JSON_STRING);
  ret->string_ = s;
  return ret;
}

JsonNode *json_mkstring(const char *s)
{
  return mkstring(json_strdup(s));
}

JsonNode *json_mknumber(double n)
{
  JsonNode *node = mknode(JSON_NUMBER);
  node->number_ = n;
  return node;
}

JsonNode *json_mkarray(void)
{
  return mknode(JSON_ARRAY);
}

JsonNode *json_mkobject(void)
{
  return mknode(JSON_OBJECT);
}

static void append_node(JsonNode *parent, JsonNode *child)
{
  if (child != NULL && parent != NULL) {
      child->parent = parent;
      child->prev = parent->children.tail;
      child->next = NULL;

      if (parent->children.tail != NULL)
          parent->children.tail->next = child;
      else
          parent->children.head = child;
      parent->children.tail = child;
  }
}

static void prepend_node(JsonNode *parent, JsonNode *child)
{
  if (child != NULL && parent != NULL) {
      child->parent = parent;
      child->prev = NULL;
      child->next = parent->children.head;

      if (parent->children.head != NULL)
          parent->children.head->prev = child;
      else
          parent->children.tail = child;
      parent->children.head = child;
  }
}

static void append_member(JsonNode *object, char *key, JsonNode *value)
{
  if (value != NULL && object != NULL) {
      value->key = key;
      append_node(object, value);
  }
}

void json_append_element(JsonNode *array, JsonNode *element)
{
  if (array != NULL && element !=NULL) {
      assert(array->tag == JSON_ARRAY);
      assert(element->parent == NULL);

      append_node(array, element);
  }
}

void json_prepend_element(JsonNode *array, JsonNode *element)
{
  assert(array->tag == JSON_ARRAY);
  assert(element->parent == NULL);

  prepend_node(array, element);
}

void json_append_member(JsonNode *object, const char *key, JsonNode *value)
{
  if (object != NULL && key != NULL && value != NULL) {
      assert(object->tag == JSON_OBJECT);
      assert(value->parent == NULL);

      append_member(object, json_strdup(key), value);
  }
}

void json_prepend_member(JsonNode *object, const char *key, JsonNode *value)
{
  if (object != NULL && key != NULL && value != NULL) {
      assert(object->tag == JSON_OBJECT);
      assert(value->parent == NULL);

      value->key = json_strdup(key);
      prepend_node(object, value);
  }
}

void json_remove_from_parent(JsonNode *node)
{
  if (node != NULL) {
      JsonNode *parent = node->parent;

      if (parent != NULL) {
          if (node->prev != NULL)
              node->prev->next = node->next;
          else
              parent->children.head = node->next;

          if (node->next != NULL)
              node->next->prev = node->prev;
          else
              parent->children.tail = node->prev;

          free(node->key);

          node->parent = NULL;
          node->prev = node->next = NULL;
          node->key = NULL;
      }
  }
}

static bool parse_value(const char **sp, JsonNode **out)
{
  const char *s = *sp;

  switch (*s) {
    case 'n':
      if (expect_literal(&s, "null")) {
        if (out)
          *out = json_mknull();
        *sp = s;
        return true;
      }
      return false;

    case 'f':
      if (expect_literal(&s, "false")) {
        if (out)
          *out = json_mkbool(false);
        *sp = s;
        return true;
      }
      return false;

    case 't':
      if (expect_literal(&s, "true")) {
        if (out)
          *out = json_mkbool(true);
        *sp = s;
        return true;
      }
      return false;

    case '"': {
      char *str = NULL;
      if (parse_string(&s, out ? &str : NULL)) {
        if (out)
          *out = mkstring(str);
        *sp = s;
        return true;
      }
      return false;
    }

    case '[':
      if (parse_array(&s, out)) {
        *sp = s;
        return true;
      }
      return false;

    case '{':
      if (parse_object(&s, out)) {
        *sp = s;
        return true;
      }
      return false;

    default: {
      double num;
      if (parse_number(&s, out ? &num : NULL)) {
        if (out)
          *out = json_mknumber(num);
        *sp = s;
        return true;
      }
      return false;
    }
  }
}

static bool parse_array(const char **sp, JsonNode **out)
{
  const char *s = *sp;
  JsonNode *ret = out ? json_mkarray() : NULL;
  JsonNode *element = NULL;

  if (*s++ != '[')
    goto failure;
  skip_space(&s);

  if (*s == ']') {
    s++;
    goto success;
  }

  for (;;) {
    if (!parse_value(&s, out ? &element : NULL))
      goto failure;
    skip_space(&s);

    if (out)
      json_append_element(ret, element);

    if (*s == ']') {
      s++;
      goto success;
    }

    if (*s++ != ',')
      goto failure;
    skip_space(&s);
  }

success:
  *sp = s;
  if (out)
    *out = ret;
  return true;

failure:
  json_delete(ret);
  return false;
}

static bool parse_object(const char **sp, JsonNode **out)
{
  const char *s = *sp;
  JsonNode *ret = out ? json_mkobject() : NULL;
  char *key = NULL;
  JsonNode *value = NULL;

  if (*s++ != '{')
    goto failure;
  skip_space(&s);

  if (*s == '}') {
    s++;
    goto success;
  }

  for (;;) {
    if (!parse_string(&s, out ? &key : NULL))
      goto failure;
    skip_space(&s);

    if (*s++ != ':')
      goto failure_free_key;
    skip_space(&s);

    if (!parse_value(&s, out ? &value : NULL))
      goto failure_free_key;
    skip_space(&s);

    if (out)
      append_member(ret, key, value);

    if (*s == '}') {
      s++;
      goto success;
    }

    if (*s++ != ',')
      goto failure;
    skip_space(&s);
  }

success:
  *sp = s;
  if (out)
    *out = ret;
  return true;

failure_free_key:
  if (out)
    free(key);
failure:
  json_delete(ret);
  return false;
}

bool parse_string(const char **sp, char **out)
{
  const char *s = *sp;
  SB sb = { 0, 0, 0 };
  char throwaway_buffer[4];
    /* enough space for a UTF-8 character */
  char *b;

  if (*s++ != '"')
    return false;

  if (out) {
    sb_init(&sb);
    sb_need(&sb, 4);
    b = sb.cur;
  } else {
    b = throwaway_buffer;
  }

  while (*s != '"') {
    unsigned char c = *s++;

    /* Parse next character, and write it to b. */
    if (c == '\\') {
      c = *s++;
      switch (c) {
        case '"':
        case '\\':
        case '/':
          *b++ = c;
          break;
        case 'b':
          *b++ = '\b';
          break;
        case 'f':
          *b++ = '\f';
          break;
        case 'n':
          *b++ = '\n';
          break;
        case 'r':
          *b++ = '\r';
          break;
        case 't':
          *b++ = '\t';
          break;
        case 'u':
        {
          uint16_t uc, lc;
          uint32_t unicode;

          if (!parse_hex16(&s, &uc))
            goto failed;

          if (uc >= 0xD800 && uc <= 0xDFFF) {
            /* Handle UTF-16 surrogate pair. */
            if (*s++ != '\\' || *s++ != 'u' || !parse_hex16(&s, &lc))
              goto failed; /* Incomplete surrogate pair. */
            if (!from_surrogate_pair(uc, lc, &unicode))
              goto failed; /* Invalid surrogate pair. */
          } else if (uc == 0) {
            /* Disallow "\u0000". */
            goto failed;
          } else {
            unicode = uc;
          }

          b += utf8_write_char(unicode, b);
          break;
        }
        default:
          /* Invalid escape */
          goto failed;
      }
    } else if (c <= 0x1F) {
      /* Control characters are not allowed in string literals. */
      goto failed;
    } else {
      /* Validate and echo a UTF-8 character. */
      int len;

      s--;
      len = utf8_validate_cz(s);
      if (len == 0)
        goto failed; /* Invalid UTF-8 character. */

      while (len--)
        *b++ = *s++;
    }

    /*
     * Update sb to know about the new bytes,
     * and set up b to write another character.
     */
    if (out) {
      sb.cur = b;
      sb_need(&sb, 4);
      b = sb.cur;
    } else {
      b = throwaway_buffer;
    }
  }
  s++;

  if (out)
    *out = sb_finish(&sb);
  *sp = s;
  return true;

failed:
  if (out)
    sb_free(&sb);
  return false;
}

bool is_space(const char *c) {
  return ((*c) == '\t' || (*c) == '\n' || (*c) == '\r' || (*c) == ' ');
}

bool is_digit(const char *c){
  return ((*c) >= '0' && (*c) <= '9');
}

/*
 * The JSON spec says that a number shall follow this precise pattern
 * (spaces and quotes added for readability):
 *   '-'? (0 | [1-9][0-9]*) ('.' [0-9]+)? ([Ee] [+-]? [0-9]+)?
 *
 * However, some JSON parsers are more liberal.  For instance, PHP accepts
 * '.5' and '1.'.  JSON.parse accepts '+3'.
 *
 * This function takes the strict approach.
 */
bool parse_number(const char **sp, double *out)
{
  const char *s = *sp;

  /* '-'? */
  if (*s == '-')
    s++;

  /* (0 | [1-9][0-9]*) */
  if (*s == '0') {
    s++;
  } else {
    if (!is_digit(s))
      return false;
    do {
      s++;
    } while (is_digit(s));
  }

  /* ('.' [0-9]+)? */
  if (*s == '.') {
    s++;
    if (!is_digit(s))
      return false;
    do {
      s++;
    } while (is_digit(s));
  }

  /* ([Ee] [+-]? [0-9]+)? */
  if (*s == 'E' || *s == 'e') {
    s++;
    if (*s == '+' || *s == '-')
      s++;
    if (!is_digit(s))
      return false;
    do {
      s++;
    } while (is_digit(s));
  }

  if (out)
    *out = strtod(*sp, NULL);

  *sp = s;
  return true;
}

static void skip_space(const char **sp)
{
  const char *s = *sp;
  while (is_space(s))
    s++;
  *sp = s;
}

static void emit_value(SB *out, const JsonNode *node)
{
  assert(tag_is_valid(node->tag));
  switch (node->tag) {
    case JSON_NULL:
      sb_puts(out, "null");
      break;
    case JSON_BOOL:
      sb_puts(out, node->bool_ ? "true" : "false");
      break;
    case JSON_STRING:
      emit_string(out, node->string_);
      break;
    case JSON_NUMBER:
      emit_number(out, node->number_);
      break;
    case JSON_ARRAY:
      emit_array(out, node);
      break;
    case JSON_OBJECT:
      emit_object(out, node);
      break;
    default:
      assert(false);
  }
}

void emit_value_indented(SB *out, const JsonNode *node, const char *space, int indent_level)
{
  assert(tag_is_valid(node->tag));
  switch (node->tag) {
    case JSON_NULL:
      sb_puts(out, "null");
      break;
    case JSON_BOOL:
      sb_puts(out, node->bool_ ? "true" : "false");
      break;
    case JSON_STRING:
      emit_string(out, node->string_);
      break;
    case JSON_NUMBER:
      emit_number(out, node->number_);
      break;
    case JSON_ARRAY:
      emit_array_indented(out, node, space, indent_level);
      break;
    case JSON_OBJECT:
      emit_object_indented(out, node, space, indent_level);
      break;
    default:
      assert(false);
  }
}

static void emit_array(SB *out, const JsonNode *array)
{
  const JsonNode *element;

  sb_putc(out, '[');
  json_foreach(element, array) {
    emit_value(out, element);
    if (element->next != NULL)
      sb_putc(out, ',');
  }
  sb_putc(out, ']');
}

static void emit_array_indented(SB *out, const JsonNode *array, const char *space, int indent_level)
{
  const JsonNode *element = array->children.head;
  int i;

  if (element == NULL) {
    sb_puts(out, "[]");
    return;
  }

  sb_puts(out, "[\n");
  while (element != NULL) {
    for (i = 0; i < indent_level + 1; i++)
      sb_puts(out, space);
    emit_value_indented(out, element, space, indent_level + 1);

    element = element->next;
    sb_puts(out, element != NULL ? ",\n" : "\n");
  }
  for (i = 0; i < indent_level; i++)
    sb_puts(out, space);
  sb_putc(out, ']');
}

static void emit_object(SB *out, const JsonNode *object)
{
  const JsonNode *member;

  sb_putc(out, '{');
  json_foreach(member, object) {
    emit_string(out, member->key);
    sb_putc(out, ':');
    emit_value(out, member);
    if (member->next != NULL)
      sb_putc(out, ',');
  }
  sb_putc(out, '}');
}

static void emit_object_indented(SB *out, const JsonNode *object, const char *space, int indent_level)
{
  const JsonNode *member = object->children.head;
  int i;

  if (member == NULL) {
    sb_puts(out, "{}");
    return;
  }

  sb_puts(out, "{\n");
  while (member != NULL) {
    for (i = 0; i < indent_level + 1; i++)
      sb_puts(out, space);
    emit_string(out, member->key);
    sb_puts(out, ": ");
    emit_value_indented(out, member, space, indent_level + 1);

    member = member->next;
    sb_puts(out, member != NULL ? ",\n" : "\n");
  }
  for (i = 0; i < indent_level; i++)
    sb_puts(out, space);
  sb_putc(out, '}');
}

void emit_string(SB *out, const char *str)
{
  bool escape_unicode = false;
  const char *s = str;
  char *b;

// make assertion catchable
#ifndef NDEBUG
  if (!utf8_validate(str)) {
    throw utf8::invalid_utf8(0);
  }
#endif

  assert(utf8_validate(str));

  /*
   * 14 bytes is enough space to write up to two
   * \uXXXX escapes and two quotation marks.
   */
  sb_need(out, 14);
  b = out->cur;

  *b++ = '"';
  while (*s != 0) {
    unsigned char c = *s++;

    /* Encode the next character, and write it to b. */
    switch (c) {
      case '"':
        *b++ = '\\';
        *b++ = '"';
        break;
      case '\\':
        *b++ = '\\';
        *b++ = '\\';
        break;
      case '\b':
        *b++ = '\\';
        *b++ = 'b';
        break;
      case '\f':
        *b++ = '\\';
        *b++ = 'f';
        break;
      case '\n':
        *b++ = '\\';
        *b++ = 'n';
        break;
      case '\r':
        *b++ = '\\';
        *b++ = 'r';
        break;
      case '\t':
        *b++ = '\\';
        *b++ = 't';
        break;
      default: {
        int len;

        s--;
        len = utf8_validate_cz(s);

        if (len == 0) {
          /*
           * Handle invalid UTF-8 character gracefully in production
           * by writing a replacement character (U+FFFD)
           * and skipping a single byte.
           *
           * This should never happen when assertions are enabled
           * due to the assertion at the beginning of this function.
           */
          assert(false);
          if (escape_unicode) {
            strcpy(b, "\\uFFFD");
            b += 6;
          } else {
            *b++ = 0xEFu;
            *b++ = 0xBFu;
            *b++ = 0xBDu;
          }
          s++;
        } else if (c < 0x1F || (c >= 0x80 && escape_unicode)) {
          /* Encode using \u.... */
          uint32_t unicode;

          s += utf8_read_char(s, &unicode);

          if (unicode <= 0xFFFF) {
            *b++ = '\\';
            *b++ = 'u';
            b += write_hex16(b, unicode);
          } else {
            /* Produce a surrogate pair. */
            uint16_t uc, lc;
            assert(unicode <= 0x10FFFF);
            to_surrogate_pair(unicode, &uc, &lc);
            *b++ = '\\';
            *b++ = 'u';
            b += write_hex16(b, uc);
            *b++ = '\\';
            *b++ = 'u';
            b += write_hex16(b, lc);
          }
        } else {
          /* Write the character directly. */
          while (len--)
            *b++ = *s++;
        }

        break;
      }
    }

    /*
     * Update *out to know about the new bytes,
     * and set up b to write another encoded character.
     */
    out->cur = b;
    sb_need(out, 14);
    b = out->cur;
  }
  *b++ = '"';

  out->cur = b;
}

static void emit_number(SB *out, double num)
{
  /*
   * This isn't exactly how JavaScript renders numbers,
   * but it should produce valid JSON for reasonable numbers
   * preserve precision well enough, and avoid some oddities
   * like 0.3 -> 0.299999999999999988898 .
   */
  char buf[64];
  sprintf(buf, "%.16g", num);

  if (number_is_valid(buf))
    sb_puts(out, buf);
  else
    sb_puts(out, "null");
}

static bool tag_is_valid(unsigned int tag)
{
  return (/* tag >= JSON_NULL && */ tag <= JSON_OBJECT);
}

static bool number_is_valid(const char *num)
{
  return (parse_number(&num, NULL) && *num == '\0');
}

static bool expect_literal(const char **sp, const char *str)
{
  const char *s = *sp;

  while (*str != '\0')
    if (*s++ != *str++)
      return false;

  *sp = s;
  return true;
}

/*
 * Parses exactly 4 hex characters (capital or lowercase).
 * Fails if any input chars are not [0-9A-Fa-f].
 */
static bool parse_hex16(const char **sp, uint16_t *out)
{
  const char *s = *sp;
  uint16_t ret = 0;
  uint16_t i;
  uint16_t tmp;
  char c;

  for (i = 0; i < 4; i++) {
    c = *s++;
    if (c >= '0' && c <= '9')
      tmp = c - '0';
    else if (c >= 'A' && c <= 'F')
      tmp = c - 'A' + 10;
    else if (c >= 'a' && c <= 'f')
      tmp = c - 'a' + 10;
    else
      return false;

    ret <<= 4;
    ret += tmp;
  }

  if (out)
    *out = ret;
  *sp = s;
  return true;
}

/*
 * Encodes a 16-bit number into hexadecimal,
 * writing exactly 4 hex chars.
 */
static int write_hex16(char *out, uint16_t val)
{
  const char *hex = "0123456789ABCDEF";

  *out++ = hex[(val >> 12) & 0xF];
  *out++ = hex[(val >> 8)  & 0xF];
  *out++ = hex[(val >> 4)  & 0xF];
  *out++ = hex[ val        & 0xF];

  return 4;
}

bool json_check(const JsonNode *node, char errmsg[256])
{
  #define problem(...) do { \
      if (errmsg != NULL) \
        snprintf(errmsg, 256, __VA_ARGS__); \
      return false; \
    } while (0)

  if (node->key != NULL && !utf8_validate(node->key))
    problem("key contains invalid UTF-8");

  if (!tag_is_valid(node->tag))
    problem("tag is invalid (%u)", node->tag);

  if (node->tag == JSON_BOOL) {
    if (node->bool_ != false && node->bool_ != true)
      problem("bool_ is neither false (%d) nor true (%d)", (int)false, (int)true);
  } else if (node->tag == JSON_STRING) {
    if (node->string_ == NULL)
      problem("string_ is NULL");
    if (!utf8_validate(node->string_))
      problem("string_ contains invalid UTF-8");
  } else if (node->tag == JSON_ARRAY || node->tag == JSON_OBJECT) {
    JsonNode *head = node->children.head;
    JsonNode *tail = node->children.tail;

    if (head == NULL || tail == NULL) {
      if (head != NULL)
        problem("tail is NULL, but head is not");
      if (tail != NULL)
        problem("head is NULL, but tail is not");
    } else {
      JsonNode *child;
      JsonNode *last = NULL;

      if (head->prev != NULL)
        problem("First child's prev pointer is not NULL");

      for (child = head; child != NULL; last = child, child = child->next) {
        if (child == node)
          problem("node is its own child");
        if (child->next == child)
          problem("child->next == child (cycle)");
        if (child->next == head)
          problem("child->next == head (cycle)");

        if (child->parent != node)
          problem("child does not point back to parent");
        if (child->next != NULL && child->next->prev != child)
          problem("child->next does not point back to child");

        if (node->tag == JSON_ARRAY && child->key != NULL)
          problem("Array element's key is not NULL");
        if (node->tag == JSON_OBJECT && child->key == NULL)
          problem("Object member's key is NULL");

        if (!json_check(child, errmsg))
          return false;
      }

      if (last != tail)
        problem("tail does not match pointer found by starting at head and following next links");
    }
  }

  return true;

  #undef problem
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.

#include <functional>
#include <locale>


namespace Sass {

  /*#########################################################################*/
  /*#########################################################################*/

  bool Selector_List::operator== (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) { return *this == *sl; }
    else if (auto ss = Cast<Simple_Selector>(&rhs)) { return *this == *ss; }
    else if (auto cpx = Cast<Complex_Selector>(&rhs)) { return *this == *cpx; }
    else if (auto cpd = Cast<Compound_Selector>(&rhs)) { return *this == *cpd; }
    else if (auto ls = Cast<List>(&rhs)) { return *this == *ls; }
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Selector_List::operator< (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) { return *this < *sl; }
    else if (auto ss = Cast<Simple_Selector>(&rhs)) { return *this < *ss; }
    else if (auto cpx = Cast<Complex_Selector>(&rhs)) { return *this < *cpx; }
    else if (auto cpd = Cast<Compound_Selector>(&rhs)) { return *this < *cpd; }
    else if (auto ls = Cast<List>(&rhs)) { return *this < *ls; }
    throw std::runtime_error("invalid selector base classes to compare");
  }

  // Selector lists can be compared to comma lists
  bool Selector_List::operator== (const Expression& rhs) const
  {
    if (auto l = Cast<List>(&rhs)) { return *this == *l; }
    if (auto s = Cast<Selector>(&rhs)) { return *this == *s; }
    if (Cast<String>(&rhs) || Cast<Null>(&rhs)) { return false; }
    throw std::runtime_error("invalid selector base classes to compare");
  }

  // Selector lists can be compared to comma lists
  bool Selector_List::operator< (const Expression& rhs) const
  {
    if (auto l = Cast<List>(&rhs)) { return *this < *l; }
    if (auto s = Cast<Selector>(&rhs)) { return *this < *s; }
    if (Cast<String>(&rhs) || Cast<Null>(&rhs)) { return true; }
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Complex_Selector::operator== (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this == *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this == *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this == *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this == *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Complex_Selector::operator< (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this < *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this < *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this < *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this < *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Compound_Selector::operator== (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this == *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this == *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this == *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this == *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Compound_Selector::operator< (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this < *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this < *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this < *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this < *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Selector_Schema::operator== (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this == *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this == *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this == *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this == *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Selector_Schema::operator< (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this < *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this < *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this < *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this < *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Simple_Selector::operator== (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this == *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this == *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this == *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this == *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  bool Simple_Selector::operator< (const Selector& rhs) const
  {
    if (auto sl = Cast<Selector_List>(&rhs)) return *this < *sl;
    if (auto ss = Cast<Simple_Selector>(&rhs)) return *this < *ss;
    if (auto cs = Cast<Complex_Selector>(&rhs)) return *this < *cs;
    if (auto ch = Cast<Compound_Selector>(&rhs)) return *this < *ch;
    throw std::runtime_error("invalid selector base classes to compare");
  }

  /*#########################################################################*/
  /*#########################################################################*/

  bool Selector_List::operator< (const Selector_List& rhs) const
  {
    size_t l = rhs.length();
    if (length() < l) l = length();
    for (size_t i = 0; i < l; i ++) {
      if (*at(i) < *rhs.at(i)) return true;
    }
    return false;
  }

  bool Selector_List::operator== (const Selector_List& rhs) const
  {
    if (&rhs == this) return true;
    if (rhs.length() != length()) return false;
    std::unordered_set<const Complex_Selector *, HashPtr, ComparePtrs> lhs_set;
    lhs_set.reserve(length());
    for (const Complex_Selector_Obj &element : elements()) {
      lhs_set.insert(element.ptr());
    }
    for (const Complex_Selector_Obj &element : rhs.elements()) {
        if (lhs_set.find(element.ptr()) == lhs_set.end()) return false;
    }
    return true;
  }

  bool Compound_Selector::operator< (const Compound_Selector& rhs) const
  {
    size_t L = std::min(length(), rhs.length());
    for (size_t i = 0; i < L; ++i)
    {
      Simple_Selector_Ptr l = (*this)[i];
      Simple_Selector_Ptr r = rhs[i];
      if (!l && !r) return false;
      else if (!r) return false;
      else if (!l) return true;
      else if (*l != *r)
      { return *l < *r; }
    }
    // just compare the length now
    return length() < rhs.length();
  }

  bool Compound_Selector::operator== (const Compound_Selector& rhs) const
  {
    if (&rhs == this) return true;
    if (rhs.length() != length()) return false;
    std::unordered_set<const Simple_Selector *, HashPtr, ComparePtrs> lhs_set;
    lhs_set.reserve(length());
    for (const Simple_Selector_Obj &element : elements()) {
      lhs_set.insert(element.ptr());
    }
    // there is no break?!
    for (const Simple_Selector_Obj &element : rhs.elements()) {
      if (lhs_set.find(element.ptr()) == lhs_set.end()) return false;
    }
    return true;
  }

  bool Complex_Selector::operator< (const Complex_Selector& rhs) const
  {
    // const iterators for tails
    Complex_Selector_Ptr_Const l = this;
    Complex_Selector_Ptr_Const r = &rhs;
    Compound_Selector_Ptr l_h = NULL;
    Compound_Selector_Ptr r_h = NULL;
    if (l) l_h = l->head();
    if (r) r_h = r->head();
    // process all tails
    while (true)
    {
      #ifdef DEBUG
      // skip empty ancestor first
      if (l && l->is_empty_ancestor())
      {
        l_h = NULL;
        l = l->tail();
        if(l) l_h = l->head();
        continue;
      }
      // skip empty ancestor first
      if (r && r->is_empty_ancestor())
      {
        r_h = NULL;
        r = r->tail();
        if (r) r_h = r->head();
        continue;
      }
      #endif
      // check for valid selectors
      if (!l) return !!r;
      if (!r) return false;
      // both are null
      else if (!l_h && !r_h)
      {
        // check combinator after heads
        if (l->combinator() != r->combinator())
        { return l->combinator() < r->combinator(); }
        // advance to next tails
        l = l->tail();
        r = r->tail();
        // fetch the next headers
        l_h = NULL; r_h = NULL;
        if (l) l_h = l->head();
        if (r) r_h = r->head();
      }
      // one side is null
      else if (!r_h) return true;
      else if (!l_h) return false;
      // heads ok and equal
      else if (*l_h == *r_h)
      {
        // check combinator after heads
        if (l->combinator() != r->combinator())
        { return l->combinator() < r->combinator(); }
        // advance to next tails
        l = l->tail();
        r = r->tail();
        // fetch the next headers
        l_h = NULL; r_h = NULL;
        if (l) l_h = l->head();
        if (r) r_h = r->head();
      }
      // heads are not equal
      else return *l_h < *r_h;
    }
  }

  bool Complex_Selector::operator== (const Complex_Selector& rhs) const
  {
    // const iterators for tails
    Complex_Selector_Ptr_Const l = this;
    Complex_Selector_Ptr_Const r = &rhs;
    Compound_Selector_Ptr l_h = NULL;
    Compound_Selector_Ptr r_h = NULL;
    if (l) l_h = l->head();
    if (r) r_h = r->head();
    // process all tails
    while (true)
    {
      // check the pointers
      if (!r) return !l;
      if (!l) return !r;
      // both are null
      if (!l_h && !r_h)
      {
        // check combinator after heads
        if (l->combinator() != r->combinator())
        { return l->combinator() < r->combinator(); }
        // advance to next tails
        l = l->tail();
        r = r->tail();
        // fetch the next heads
        l_h = NULL; r_h = NULL;
        if (l) l_h = l->head();
        if (r) r_h = r->head();
      }
      // equals if other head is empty
      else if ((!l_h && !r_h) ||
               (!l_h && r_h->empty()) ||
               (!r_h && l_h->empty()) ||
               (l_h && r_h && *l_h == *r_h))
      {
        // check combinator after heads
        if (l->combinator() != r->combinator())
        { return l->combinator() == r->combinator(); }
        // advance to next tails
        l = l->tail();
        r = r->tail();
        // fetch the next heads
        l_h = NULL; r_h = NULL;
        if (l) l_h = l->head();
        if (r) r_h = r->head();
      }
      // abort
      else break;
    }
    // unreachable
    return false;
  }

  /*#########################################################################*/
  /*#########################################################################*/

  bool Selector_List::operator== (const Complex_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return rhs.empty();
    return *at(0) == rhs;
  }

  bool Selector_List::operator< (const Complex_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return !rhs.empty();
    return *at(0) < rhs;
  }

  bool Selector_List::operator== (const Compound_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return rhs.empty();
    return *at(0) == rhs;
  }

  bool Selector_List::operator< (const Compound_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return !rhs.empty();
    return *at(0) < rhs;
  }

  bool Selector_List::operator== (const Simple_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return rhs.empty();
    return *at(0) == rhs;
  }

  bool Selector_List::operator< (const Simple_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return !rhs.empty();
    return *at(0) < rhs;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Complex_Selector::operator== (const Selector_List& rhs) const
  {
    size_t len = rhs.length();
    if (len > 1) return false;
    if (len == 0) return empty();
    return *this == *rhs.at(0);
  }

  bool Complex_Selector::operator< (const Selector_List& rhs) const
  {
    size_t len = rhs.length();
    if (len > 1) return true;
    if (len == 0) return false;
    return *this < *rhs.at(0);
  }

  bool Complex_Selector::operator== (const Compound_Selector& rhs) const
  {
    if (tail()) return false;
    if (!head()) return rhs.empty();
    return *head() == rhs;
  }

  bool Complex_Selector::operator< (const Compound_Selector& rhs) const
  {
    if (tail()) return false;
    if (!head()) return !rhs.empty();
    return *head() < rhs;
  }

  bool Complex_Selector::operator== (const Simple_Selector& rhs) const
  {
    if (tail()) return false;
    if (!head()) return rhs.empty();
    return *head() == rhs;
  }

  bool Complex_Selector::operator< (const Simple_Selector& rhs) const
  {
    if (tail()) return false;
    if (!head()) return !rhs.empty();
    return *head() < rhs;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Compound_Selector::operator== (const Selector_List& rhs) const
  {
    size_t len = rhs.length();
    if (len > 1) return false;
    if (len == 0) return empty();
    return *this == *rhs.at(0);
  }

  bool Compound_Selector::operator< (const Selector_List& rhs) const
  {
    size_t len = rhs.length();
    if (len > 1) return true;
    if (len == 0) return false;
    return *this < *rhs.at(0);
  }

  bool Compound_Selector::operator== (const Complex_Selector& rhs) const
  {
    if (rhs.tail()) return false;
    if (!rhs.head()) return empty();
    return *this == *rhs.head();
  }

  bool Compound_Selector::operator< (const Complex_Selector& rhs) const
  {
    if (rhs.tail()) return true;
    if (!rhs.head()) return false;
    return *this < *rhs.head();
  }

  bool Compound_Selector::operator< (const Simple_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return rhs.empty();
    return *at(0) == rhs;
  }

  bool Compound_Selector::operator== (const Simple_Selector& rhs) const
  {
    size_t len = length();
    if (len > 1) return false;
    if (len == 0) return !rhs.empty();
    return *at(0) < rhs;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Simple_Selector::operator== (const Selector_List& rhs) const
  {
    return rhs.length() == 1 && *this == *rhs.at(0);
  }

  bool Simple_Selector::operator< (const Selector_List& rhs) const
  {
    size_t len = rhs.length();
    if (len > 1) return true;
    if (len == 0) return false;
    return *this < *rhs.at(0);
  }

  bool Simple_Selector::operator== (const Complex_Selector& rhs) const
  {
    return !rhs.tail() && rhs.head() &&
           rhs.combinator() == Complex_Selector::ANCESTOR_OF &&
           *this == *rhs.head();
  }

  bool Simple_Selector::operator< (const Complex_Selector& rhs) const
  {
    if (rhs.tail()) return true;
    if (!rhs.head()) return false;
    return *this < *rhs.head();
  }

  bool Simple_Selector::operator== (const Compound_Selector& rhs) const
  {
    return rhs.length() == 1 && *this == *rhs.at(0);
  }

  bool Simple_Selector::operator< (const Compound_Selector& rhs) const
  {
    size_t len = rhs.length();
    if (len > 1) return true;
    if (len == 0) return false;
    return *this < *rhs.at(0);
  }

  /*#########################################################################*/
  /*#########################################################################*/

  bool Simple_Selector::operator== (const Simple_Selector& rhs) const
  {
    switch (simple_type()) {
      case ID_SEL: return (const Id_Selector&) *this == rhs; break;
      case TYPE_SEL: return (const Type_Selector&) *this == rhs; break;
      case CLASS_SEL: return (const Class_Selector&) *this == rhs; break;
      case PARENT_SEL: return (const Parent_Selector&) *this == rhs; break;
      case PSEUDO_SEL: return (const Pseudo_Selector&) *this == rhs; break;
      case WRAPPED_SEL: return (const Wrapped_Selector&) *this == rhs; break;
      case ATTRIBUTE_SEL: return (const Attribute_Selector&) *this == rhs; break;
      case PLACEHOLDER_SEL: return (const Placeholder_Selector&) *this == rhs; break;
    }
    return false;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Id_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Id_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Type_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Type_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Class_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Class_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Parent_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Parent_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Pseudo_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Pseudo_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Wrapped_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Wrapped_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Attribute_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Attribute_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  bool Placeholder_Selector::operator== (const Simple_Selector& rhs) const
  {
    auto sel = Cast<Placeholder_Selector>(&rhs);
    return sel ? *this == *sel : false;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Id_Selector::operator== (const Id_Selector& rhs) const
  {
    // ID has no namespacing
    return name() == rhs.name();
  }

  bool Type_Selector::operator== (const Type_Selector& rhs) const
  {
    return is_ns_eq(rhs) && name() == rhs.name();
  }

  bool Class_Selector::operator== (const Class_Selector& rhs) const
  {
    // Class has no namespacing
    return name() == rhs.name();
  }

  bool Parent_Selector::operator== (const Parent_Selector& rhs) const
  {
    // Parent has no namespacing
    return name() == rhs.name();
  }

  bool Pseudo_Selector::operator== (const Pseudo_Selector& rhs) const
  {
    std::string lname = name();
    std::string rname = rhs.name();
    if (is_pseudo_class_element(lname)) {
      if (rname[0] == ':' && rname[1] == ':') {
        lname = lname.substr(1, std::string::npos);
      }
    }
    // right hand is special pseudo (single colon)
    if (is_pseudo_class_element(rname)) {
      if (lname[0] == ':' && lname[1] == ':') {
        lname = lname.substr(1, std::string::npos);
      }
    }
    // Pseudo has no namespacing
    if (lname != rname) return false;
    String_Obj lhs_ex = expression();
    String_Obj rhs_ex = rhs.expression();
    if (rhs_ex && lhs_ex) return *lhs_ex == *rhs_ex;
    else return lhs_ex.ptr() == rhs_ex.ptr();
  }

  bool Wrapped_Selector::operator== (const Wrapped_Selector& rhs) const
  {
    // Wrapped has no namespacing
    if (name() != rhs.name()) return false;
    return *(selector()) == *(rhs.selector());
  }

  bool Attribute_Selector::operator== (const Attribute_Selector& rhs) const
  {
    // get optional value state
    bool no_lhs_val = value().isNull();
    bool no_rhs_val = rhs.value().isNull();
    // both are null, therefore equal
    if (no_lhs_val && no_rhs_val) {
      return (name() == rhs.name())
        && (matcher() == rhs.matcher())
        && (is_ns_eq(rhs));
    }
    // both are defined, evaluate
    if (no_lhs_val == no_rhs_val) {
      return (name() == rhs.name())
        && (matcher() == rhs.matcher())
        && (is_ns_eq(rhs))
        && (*value() == *rhs.value());
    }
    // not equal
    return false;
  }

  bool Placeholder_Selector::operator== (const Placeholder_Selector& rhs) const
  {
    // Placeholder has no namespacing
    return name() == rhs.name();
  }

  /*#########################################################################*/
  /*#########################################################################*/

  bool Simple_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (simple_type()) {
      case ID_SEL: return (const Id_Selector&) *this < rhs; break;
      case TYPE_SEL: return (const Type_Selector&) *this < rhs; break;
      case CLASS_SEL: return (const Class_Selector&) *this < rhs; break;
      case PARENT_SEL: return (const Parent_Selector&) *this < rhs; break;
      case PSEUDO_SEL: return (const Pseudo_Selector&) *this < rhs; break;
      case WRAPPED_SEL: return (const Wrapped_Selector&) *this < rhs; break;
      case ATTRIBUTE_SEL: return (const Attribute_Selector&) *this < rhs; break;
      case PLACEHOLDER_SEL: return (const Placeholder_Selector&) *this < rhs; break;
    }
    return false;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Id_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case TYPE_SEL: return '#' < 's'; break;
      case CLASS_SEL: return '#' < '.'; break;
      case PARENT_SEL: return '#' < '&'; break;
      case PSEUDO_SEL: return '#' < ':'; break;
      case WRAPPED_SEL: return '#' < '('; break;
      case ATTRIBUTE_SEL: return '#' < '['; break;
      case PLACEHOLDER_SEL: return '#' < '%'; break;
      case ID_SEL: /* let if fall through */ break;
    }
    const Id_Selector& sel =
      (const Id_Selector&) rhs;
    return *this < sel;
  }

  bool Type_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return 'e' < '#'; break;
      case CLASS_SEL: return 'e' < '.'; break;
      case PARENT_SEL: return 'e' < '&'; break;
      case PSEUDO_SEL: return 'e' < ':'; break;
      case WRAPPED_SEL: return 'e' < '('; break;
      case ATTRIBUTE_SEL: return 'e' < '['; break;
      case PLACEHOLDER_SEL: return 'e' < '%'; break;
      case TYPE_SEL: /* let if fall through */ break;
    }
    const Type_Selector& sel =
      (const Type_Selector&) rhs;
    return *this < sel;
  }

  bool Class_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return '.' < '#'; break;
      case TYPE_SEL: return '.' < 's'; break;
      case PARENT_SEL: return '.' < '&'; break;
      case PSEUDO_SEL: return '.' < ':'; break;
      case WRAPPED_SEL: return '.' < '('; break;
      case ATTRIBUTE_SEL: return '.' < '['; break;
      case PLACEHOLDER_SEL: return '.' < '%'; break;
      case CLASS_SEL: /* let if fall through */ break;
    }
    const Class_Selector& sel =
      (const Class_Selector&) rhs;
    return *this < sel;
  }

  bool Pseudo_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return ':' < '#'; break;
      case TYPE_SEL: return ':' < 's'; break;
      case CLASS_SEL: return ':' < '.'; break;
      case PARENT_SEL: return ':' < '&'; break;
      case WRAPPED_SEL: return ':' < '('; break;
      case ATTRIBUTE_SEL: return ':' < '['; break;
      case PLACEHOLDER_SEL: return ':' < '%'; break;
      case PSEUDO_SEL: /* let if fall through */ break;
    }
    const Pseudo_Selector& sel =
      (const Pseudo_Selector&) rhs;
    return *this < sel;
  }

  bool Wrapped_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return '(' < '#'; break;
      case TYPE_SEL: return '(' < 's'; break;
      case CLASS_SEL: return '(' < '.'; break;
      case PARENT_SEL: return '(' < '&'; break;
      case PSEUDO_SEL: return '(' < ':'; break;
      case ATTRIBUTE_SEL: return '(' < '['; break;
      case PLACEHOLDER_SEL: return '(' < '%'; break;
      case WRAPPED_SEL: /* let if fall through */ break;
    }
    const Wrapped_Selector& sel =
      (const Wrapped_Selector&) rhs;
    return *this < sel;
  }

  bool Parent_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return '&' < '#'; break;
      case TYPE_SEL: return '&' < 's'; break;
      case CLASS_SEL: return '&' < '.'; break;
      case PSEUDO_SEL: return '&' < ':'; break;
      case WRAPPED_SEL: return '&' < '('; break;
      case ATTRIBUTE_SEL: return '&' < '['; break;
      case PLACEHOLDER_SEL: return '&' < '%'; break;
      case PARENT_SEL: /* let if fall through */ break;
    }
    const Parent_Selector& sel =
      (const Parent_Selector&) rhs;
    return *this < sel;
  }

  bool Attribute_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return '[' < '#'; break;
      case TYPE_SEL: return '[' < 'e'; break;
      case CLASS_SEL: return '[' < '.'; break;
      case PARENT_SEL: return '[' < '&'; break;
      case PSEUDO_SEL: return '[' < ':'; break;
      case WRAPPED_SEL: return '[' < '('; break;
      case PLACEHOLDER_SEL: return '[' < '%'; break;
      case ATTRIBUTE_SEL: /* let if fall through */ break;
    }
    const Attribute_Selector& sel =
      (const Attribute_Selector&) rhs;
    return *this < sel;
  }

  bool Placeholder_Selector::operator< (const Simple_Selector& rhs) const
  {
    switch (rhs.simple_type()) {
      case ID_SEL: return '%' < '#'; break;
      case TYPE_SEL: return '%' < 's'; break;
      case CLASS_SEL: return '%' < '.'; break;
      case PARENT_SEL: return '%' < '&'; break;
      case PSEUDO_SEL: return '%' < ':'; break;
      case WRAPPED_SEL: return '%' < '('; break;
      case ATTRIBUTE_SEL: return '%' < '['; break;
      case PLACEHOLDER_SEL: /* let if fall through */ break;
    }
    const Placeholder_Selector& sel =
      (const Placeholder_Selector&) rhs;
    return *this < sel;
  }

  /***************************************************************************/
  /***************************************************************************/

  bool Id_Selector::operator< (const Id_Selector& rhs) const
  {
    // ID has no namespacing
    return name() < rhs.name();
  }

  bool Type_Selector::operator< (const Type_Selector& rhs) const
  {
    return has_ns_ == rhs.has_ns_
      ? (ns_ == rhs.ns_
         ? name_ < rhs.name_
         : ns_ < rhs.ns_)
      : has_ns_ < rhs.has_ns_;
  }

  bool Class_Selector::operator< (const Class_Selector& rhs) const
  {
    // Class has no namespacing
    return name() < rhs.name();
  }

  bool Parent_Selector::operator< (const Parent_Selector& rhs) const
  {
    // Parent has no namespacing
    return name() < rhs.name();
  }

  bool Pseudo_Selector::operator< (const Pseudo_Selector& rhs) const
  {
    std::string lname = name();
    std::string rname = rhs.name();
    if (is_pseudo_class_element(lname)) {
      if (rname[0] == ':' && rname[1] == ':') {
        lname = lname.substr(1, std::string::npos);
      }
    }
    // right hand is special pseudo (single colon)
    if (is_pseudo_class_element(rname)) {
      if (lname[0] == ':' && lname[1] == ':') {
        lname = lname.substr(1, std::string::npos);
      }
    }
    // Peudo has no namespacing
    if (lname != rname)
    { return lname < rname; }
    String_Obj lhs_ex = expression();
    String_Obj rhs_ex = rhs.expression();
    if (rhs_ex && lhs_ex) return *lhs_ex < *rhs_ex;
    else return lhs_ex.ptr() < rhs_ex.ptr();
  }

  bool Wrapped_Selector::operator< (const Wrapped_Selector& rhs) const
  {
    // Wrapped has no namespacing
    if (name() != rhs.name())
    { return name() < rhs.name(); }
    return *(selector()) < *(rhs.selector());
  }

  bool Attribute_Selector::operator< (const Attribute_Selector& rhs) const
  {
    if (is_ns_eq(rhs)) {
      if (name() != rhs.name())
      { return name() < rhs.name(); }
      if (matcher() != rhs.matcher())
      { return matcher() < rhs.matcher(); }
      bool no_lhs_val = value().isNull();
      bool no_rhs_val = rhs.value().isNull();
      if (no_lhs_val && no_rhs_val) return false; // equal
      else if (no_lhs_val) return true; // lhs is null
      else if (no_rhs_val) return false; // rhs is null
      return *value() < *rhs.value(); // both are given
    } else { return ns() < rhs.ns(); }
  }

  bool Placeholder_Selector::operator< (const Placeholder_Selector& rhs) const
  {
    // Placeholder has no namespacing
    return name() < rhs.name();
  }

  /*#########################################################################*/
  /*#########################################################################*/

}

namespace Sass {

  namespace Functions {

    Signature selector_nest_sig = "selector-nest($selectors...)";
    BUILT_IN(selector_nest)
    {
      List_Ptr arglist = ARG("$selectors", List);

      // Not enough parameters
      if( arglist->length() == 0 )
        error("$selectors: At least one selector must be passed for `selector-nest'", pstate, traces);

      // Parse args into vector of selectors
      SelectorStack parsedSelectors;
      for (size_t i = 0, L = arglist->length(); i < L; ++i) {
        Expression_Obj exp = Cast<Expression>(arglist->value_at_index(i));
        if (exp->concrete_type() == Expression::NULL_VAL) {
          std::stringstream msg;
          msg << "$selectors: null is not a valid selector: it must be a string,\n";
          msg << "a list of strings, or a list of lists of strings for 'selector-nest'";
          error(msg.str(), pstate, traces);
        }
        if (String_Constant_Obj str = Cast<String_Constant>(exp)) {
          str->quote_mark(0);
        }
        std::string exp_src = exp->to_string(ctx.c_options);
        Selector_List_Obj sel = Parser::parse_selector(exp_src.c_str(), ctx, traces);
        parsedSelectors.push_back(sel);
      }

      // Nothing to do
      if( parsedSelectors.empty() ) {
        return SASS_MEMORY_NEW(Null, pstate);
      }

      // Set the first element as the `result`, keep appending to as we go down the parsedSelector vector.
      SelectorStack::iterator itr = parsedSelectors.begin();
      Selector_List_Obj result = *itr;
      ++itr;

      for(;itr != parsedSelectors.end(); ++itr) {
        Selector_List_Obj child = *itr;
        std::vector<Complex_Selector_Obj> exploded;
        selector_stack.push_back(result);
        Selector_List_Obj rv = child->resolve_parent_refs(selector_stack, traces);
        selector_stack.pop_back();
        for (size_t m = 0, mLen = rv->length(); m < mLen; ++m) {
          exploded.push_back((*rv)[m]);
        }
        result->elements(exploded);
      }

      Listize listize;
      return Cast<Value>(result->perform(&listize));
    }

    Signature selector_append_sig = "selector-append($selectors...)";
    BUILT_IN(selector_append)
    {
      List_Ptr arglist = ARG("$selectors", List);

      // Not enough parameters
      if( arglist->length() == 0 )
        error("$selectors: At least one selector must be passed for `selector-append'", pstate, traces);

      // Parse args into vector of selectors
      SelectorStack parsedSelectors;
      for (size_t i = 0, L = arglist->length(); i < L; ++i) {
        Expression_Obj exp = Cast<Expression>(arglist->value_at_index(i));
        if (exp->concrete_type() == Expression::NULL_VAL) {
          std::stringstream msg;
          msg << "$selectors: null is not a valid selector: it must be a string,\n";
          msg << "a list of strings, or a list of lists of strings for 'selector-append'";
          error(msg.str(), pstate, traces);
        }
        if (String_Constant_Ptr str = Cast<String_Constant>(exp)) {
          str->quote_mark(0);
        }
        std::string exp_src = exp->to_string();
        Selector_List_Obj sel = Parser::parse_selector(exp_src.c_str(), ctx, traces,
                                                       exp->pstate(), pstate.src,
                                                       /*allow_parent=*/false);

        parsedSelectors.push_back(sel);
      }

      // Nothing to do
      if( parsedSelectors.empty() ) {
        return SASS_MEMORY_NEW(Null, pstate);
      }

      // Set the first element as the `result`, keep appending to as we go down the parsedSelector vector.
      SelectorStack::iterator itr = parsedSelectors.begin();
      Selector_List_Obj result = *itr;
      ++itr;

      for(;itr != parsedSelectors.end(); ++itr) {
        Selector_List_Obj child = *itr;
        std::vector<Complex_Selector_Obj> newElements;

        // For every COMPLEX_SELECTOR in `result`
        // For every COMPLEX_SELECTOR in `child`
          // let parentSeqClone equal a copy of result->elements[i]
          // let childSeq equal child->elements[j]
          // Append all of childSeq head elements into parentSeqClone
          // Set the innermost tail of parentSeqClone, to childSeq's tail
        // Replace result->elements with newElements
        for (size_t i = 0, resultLen = result->length(); i < resultLen; ++i) {
          for (size_t j = 0, childLen = child->length(); j < childLen; ++j) {
            Complex_Selector_Obj parentSeqClone = SASS_MEMORY_CLONE((*result)[i]);
            Complex_Selector_Obj childSeq = (*child)[j];
            Complex_Selector_Obj base = childSeq->tail();

            // Must be a simple sequence
            if( childSeq->combinator() != Complex_Selector::Combinator::ANCESTOR_OF ) {
              error("Can't append \"" + childSeq->to_string() + "\" to \"" +
                parentSeqClone->to_string() + "\" for `selector-append'", pstate, traces);
            }

            // Cannot be a Universal selector
            Type_Selector_Obj pType = Cast<Type_Selector>(childSeq->head()->first());
            if(pType && pType->name() == "*") {
              error("Can't append \"" + childSeq->to_string() + "\" to \"" +
                parentSeqClone->to_string() + "\" for `selector-append'", pstate, traces);
            }

            // TODO: Add check for namespace stuff

            Complex_Selector_Ptr lastComponent = parentSeqClone->mutable_last();
            if (lastComponent->head() == nullptr) {
              std::string msg = "Parent \"" + parentSeqClone->to_string() + "\" is incompatible with \"" + base->to_string() + "\"";
              error(msg, pstate, traces);
            }
            lastComponent->head()->concat(base->head());
            lastComponent->tail(base->tail());

            newElements.push_back(parentSeqClone);
          }
        }

        result->elements(newElements);
      }

      Listize listize;
      return Cast<Value>(result->perform(&listize));
    }

    Signature selector_unify_sig = "selector-unify($selector1, $selector2)";
    BUILT_IN(selector_unify)
    {
      Selector_List_Obj selector1 = ARGSELS("$selector1");
      Selector_List_Obj selector2 = ARGSELS("$selector2");

      Selector_List_Obj result = selector1->unify_with(selector2);
      Listize listize;
      return Cast<Value>(result->perform(&listize));
    }

    Signature simple_selectors_sig = "simple-selectors($selector)";
    BUILT_IN(simple_selectors)
    {
      Compound_Selector_Obj sel = ARGSEL("$selector");

      List_Ptr l = SASS_MEMORY_NEW(List, sel->pstate(), sel->length(), SASS_COMMA);

      for (size_t i = 0, L = sel->length(); i < L; ++i) {
        Simple_Selector_Obj ss = (*sel)[i];
        std::string ss_string = ss->to_string() ;

        l->append(SASS_MEMORY_NEW(String_Quoted, ss->pstate(), ss_string));
      }

      return l;
    }

    Signature selector_extend_sig = "selector-extend($selector, $extendee, $extender)";
    BUILT_IN(selector_extend)
    {
      Selector_List_Obj  selector = ARGSELS("$selector");
      Selector_List_Obj  extendee = ARGSELS("$extendee");
      Selector_List_Obj  extender = ARGSELS("$extender");

      Subset_Map subset_map;
      extender->populate_extends(extendee, subset_map);
      Extend extend(subset_map);

      Selector_List_Obj result = extend.extendSelectorList(selector, false);

      Listize listize;
      return Cast<Value>(result->perform(&listize));
    }

    Signature selector_replace_sig = "selector-replace($selector, $original, $replacement)";
    BUILT_IN(selector_replace)
    {
      Selector_List_Obj selector = ARGSELS("$selector");
      Selector_List_Obj original = ARGSELS("$original");
      Selector_List_Obj replacement = ARGSELS("$replacement");
      Subset_Map subset_map;
      replacement->populate_extends(original, subset_map);
      Extend extend(subset_map);

      Selector_List_Obj result = extend.extendSelectorList(selector, true);

      Listize listize;
      return Cast<Value>(result->perform(&listize));
    }

    Signature selector_parse_sig = "selector-parse($selector)";
    BUILT_IN(selector_parse)
    {
      Selector_List_Obj sel = ARGSELS("$selector");

      Listize listize;
      return Cast<Value>(sel->perform(&listize));
    }

    Signature is_superselector_sig = "is-superselector($super, $sub)";
    BUILT_IN(is_superselector)
    {
      Selector_List_Obj  sel_sup = ARGSELS("$super");
      Selector_List_Obj  sel_sub = ARGSELS("$sub");
      bool result = sel_sup->is_superselector_of(sel_sub);
      return SASS_MEMORY_NEW(Boolean, pstate, result);
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  // Custom_Error is a valid value
  Value_Ptr To_Value::operator()(Custom_Error_Ptr e)
  {
    return e;
  }

  // Custom_Warning is a valid value
  Value_Ptr To_Value::operator()(Custom_Warning_Ptr w)
  {
    return w;
  }

  // Boolean is a valid value
  Value_Ptr To_Value::operator()(Boolean_Ptr b)
  {
    return b;
  }

  // Number is a valid value
  Value_Ptr To_Value::operator()(Number_Ptr n)
  {
    return n;
  }

  // Color is a valid value
  Value_Ptr To_Value::operator()(Color_RGBA_Ptr c)
  {
    return c;
  }

  // Color is a valid value
  Value_Ptr To_Value::operator()(Color_HSLA_Ptr c)
  {
    return c;
  }

  // String_Constant is a valid value
  Value_Ptr To_Value::operator()(String_Constant_Ptr s)
  {
    return s;
  }

  // String_Quoted is a valid value
  Value_Ptr To_Value::operator()(String_Quoted_Ptr s)
  {
    return s;
  }

  // List is a valid value
  Value_Ptr To_Value::operator()(List_Ptr l)
  {
    List_Obj ll = SASS_MEMORY_NEW(List,
                               l->pstate(),
                               l->length(),
                               l->separator(),
                               l->is_arglist(),
                               l->is_bracketed());
    for (size_t i = 0, L = l->length(); i < L; ++i) {
      ll->append((*l)[i]->perform(this));
    }
    return ll.detach();
  }

  // Map is a valid value
  Value_Ptr To_Value::operator()(Map_Ptr m)
  {
    return m;
  }

  // Null is a valid value
  Value_Ptr To_Value::operator()(Null_Ptr n)
  {
    return n;
  }

  // Function is a valid value
  Value_Ptr To_Value::operator()(Function_Ptr n)
  {
    return n;
  }

  // Argument returns its value
  Value_Ptr To_Value::operator()(Argument_Ptr arg)
  {
    if (!arg->name().empty()) return 0;
    return arg->value()->perform(this);
  }

  // Selector_List is converted to a string
  Value_Ptr To_Value::operator()(Selector_List_Ptr s)
  {
    return SASS_MEMORY_NEW(String_Quoted,
                           s->pstate(),
                           s->to_string(ctx.c_options));
  }

  // Binary_Expression is converted to a string
  Value_Ptr To_Value::operator()(Binary_Expression_Ptr s)
  {
    return SASS_MEMORY_NEW(String_Quoted,
                           s->pstate(),
                           s->to_string(ctx.c_options));
  }

};


#ifdef DEBUG_SHARED_PTR
#endif

namespace Sass {

  #ifdef DEBUG_SHARED_PTR
  void SharedObj::dumpMemLeaks() {
    if (!all.empty()) {
      std::cerr << "###################################\n";
      std::cerr << "# REPORTING MISSING DEALLOCATIONS #\n";
      std::cerr << "###################################\n";
      for (SharedObj* var : all) {
        if (AST_Node_Ptr ast = dynamic_cast<AST_Node*>(var)) {
          debug_ast(ast);
        } else {
          std::cerr << "LEAKED " << var << "\n";
        }
      }
    }
  }
  std::vector<SharedObj*> SharedObj::all;
  #endif

  bool SharedObj::taint = false;
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


// Notes about delayed: some ast nodes can have delayed evaluation so
// they can preserve their original semantics if needed. This is most
// prominently exhibited by the division operation, since it is not
// only a valid operation, but also a valid css statement (i.e. for
// fonts, as in `16px/24px`). When parsing lists and expression we
// unwrap single items from lists and other operations. A nested list
// must not be delayed, only the items of the first level sometimes
// are delayed (as with argument lists). To achieve this we need to
// pass status to the list parser, so this can be set correctly.
// Another case with delayed values are colors. In compressed mode
// only processed values get compressed (other are left as written).


namespace Sass {
  using namespace Constants;
  using namespace Prelexer;

  Parser Parser::from_c_str(const char* beg, Context& ctx, Backtraces traces, ParserState pstate, const char* source, bool allow_parent)
  {
    pstate.offset.column = 0;
    pstate.offset.line = 0;
    Parser p(ctx, pstate, traces, allow_parent);
    p.source   = source ? source : beg;
    p.position = beg ? beg : p.source;
    p.end      = p.position + strlen(p.position);
    Block_Obj root = SASS_MEMORY_NEW(Block, pstate);
    p.block_stack.push_back(root);
    root->is_root(true);
    return p;
  }

  Parser Parser::from_c_str(const char* beg, const char* end, Context& ctx, Backtraces traces, ParserState pstate, const char* source, bool allow_parent)
  {
    pstate.offset.column = 0;
    pstate.offset.line = 0;
    Parser p(ctx, pstate, traces, allow_parent);
    p.source   = source ? source : beg;
    p.position = beg ? beg : p.source;
    p.end      = end ? end : p.position + strlen(p.position);
    Block_Obj root = SASS_MEMORY_NEW(Block, pstate);
    p.block_stack.push_back(root);
    root->is_root(true);
    return p;
  }

   void Parser::advanceToNextToken() {
      lex < css_comments >(false);
      // advance to position
      pstate += pstate.offset;
      pstate.offset.column = 0;
      pstate.offset.line = 0;
    }

  Selector_List_Obj Parser::parse_selector(const char* beg, Context& ctx, Backtraces traces, ParserState pstate, const char* source, bool allow_parent)
  {
    Parser p = Parser::from_c_str(beg, ctx, traces, pstate, source, allow_parent);
    // ToDo: remap the source-map entries somehow
    return p.parse_selector_list(false);
  }

  bool Parser::peek_newline(const char* start)
  {
    return peek_linefeed(start ? start : position)
           && ! peek_css<exactly<'{'>>(start);
  }

  Parser Parser::from_token(Token t, Context& ctx, Backtraces traces, ParserState pstate, const char* source)
  {
    Parser p(ctx, pstate, traces);
    p.source   = source ? source : t.begin;
    p.position = t.begin ? t.begin : p.source;
    p.end      = t.end ? t.end : p.position + strlen(p.position);
    Block_Obj root = SASS_MEMORY_NEW(Block, pstate);
    p.block_stack.push_back(root);
    root->is_root(true);
    return p;
  }

  /* main entry point to parse root block */
  Block_Obj Parser::parse()
  {

    // consume unicode BOM
    read_bom();

    // scan the input to find invalid utf8 sequences
    const char* it = utf8::find_invalid(position, end);

    // report invalid utf8
    if (it != end) {
      pstate += Offset::init(position, it);
      traces.push_back(Backtrace(pstate));
      throw Exception::InvalidSass(pstate, traces, "Invalid UTF-8 sequence");
    }

    // create a block AST node to hold children
    Block_Obj root = SASS_MEMORY_NEW(Block, pstate, 0, true);

    // check seems a bit esoteric but works
    if (ctx.resources.size() == 1) {
      // apply headers only on very first include
      ctx.apply_custom_headers(root, path, pstate);
    }

    // parse children nodes
    block_stack.push_back(root);
    parse_block_nodes(true);
    block_stack.pop_back();

    // update final position
    root->update_pstate(pstate);

    if (position != end) {
      css_error("Invalid CSS", " after ", ": expected selector or at-rule, was ");
    }

    return root;
  }


  // convenience function for block parsing
  // will create a new block ad-hoc for you
  // this is the base block parsing function
  Block_Obj Parser::parse_css_block(bool is_root)
  {

    // parse comments before block
    // lex < optional_css_comments >();

    // lex mandatory opener or error out
    if (!lex_css < exactly<'{'> >()) {
      css_error("Invalid CSS", " after ", ": expected \"{\", was ");
    }
    // create new block and push to the selector stack
    Block_Obj block = SASS_MEMORY_NEW(Block, pstate, 0, is_root);
    block_stack.push_back(block);

    if (!parse_block_nodes(is_root)) css_error("Invalid CSS", " after ", ": expected \"}\", was ");

    if (!lex_css < exactly<'}'> >()) {
      css_error("Invalid CSS", " after ", ": expected \"}\", was ");
    }

    // update for end position
    // this seems to be done somewhere else
    // but that fixed selector schema issue
    // block->update_pstate(pstate);

    // parse comments after block
    // lex < optional_css_comments >();

    block_stack.pop_back();

    return block;
  }

  // convenience function for block parsing
  // will create a new block ad-hoc for you
  // also updates the `in_at_root` flag
  Block_Obj Parser::parse_block(bool is_root)
  {
    return parse_css_block(is_root);
  }

  // the main block parsing function
  // parses stuff between `{` and `}`
  bool Parser::parse_block_nodes(bool is_root)
  {

    // loop until end of string
    while (position < end) {

      // we should be able to refactor this
      parse_block_comments();
      lex < css_whitespace >();

      if (lex < exactly<';'> >()) continue;
      if (peek < end_of_file >()) return true;
      if (peek < exactly<'}'> >()) return true;

      if (parse_block_node(is_root)) continue;

      parse_block_comments();

      if (lex_css < exactly<';'> >()) continue;
      if (peek_css < end_of_file >()) return true;
      if (peek_css < exactly<'}'> >()) return true;

      // illegal sass
      return false;
    }
    // return success
    return true;
  }

  // parser for a single node in a block
  // semicolons must be lexed beforehand
  bool Parser::parse_block_node(bool is_root) {

    Block_Obj block = block_stack.back();

    parse_block_comments();

    // throw away white-space
    // includes line comments
    lex < css_whitespace >();

    Lookahead lookahead_result;

    // also parse block comments

    // first parse everything that is allowed in functions
    if (lex < variable >(true)) { block->append(parse_assignment()); }
    else if (lex < kwd_err >(true)) { block->append(parse_error()); }
    else if (lex < kwd_dbg >(true)) { block->append(parse_debug()); }
    else if (lex < kwd_warn >(true)) { block->append(parse_warning()); }
    else if (lex < kwd_if_directive >(true)) { block->append(parse_if_directive()); }
    else if (lex < kwd_for_directive >(true)) { block->append(parse_for_directive()); }
    else if (lex < kwd_each_directive >(true)) { block->append(parse_each_directive()); }
    else if (lex < kwd_while_directive >(true)) { block->append(parse_while_directive()); }
    else if (lex < kwd_return_directive >(true)) { block->append(parse_return_directive()); }

    // parse imports to process later
    else if (lex < kwd_import >(true)) {
      Scope parent = stack.empty() ? Scope::Rules : stack.back();
      if (parent != Scope::Function && parent != Scope::Root && parent != Scope::Rules && parent != Scope::Media) {
        if (! peek_css< uri_prefix >(position)) { // this seems to go in ruby sass 3.4.20
          error("Import directives may not be used within control directives or mixins.");
        }
      }
      // this puts the parsed doc into sheets
      // import stub will fetch this in expand
      Import_Obj imp = parse_import();
      // if it is a url, we only add the statement
      if (!imp->urls().empty()) block->append(imp);
      // process all resources now (add Import_Stub nodes)
      for (size_t i = 0, S = imp->incs().size(); i < S; ++i) {
        block->append(SASS_MEMORY_NEW(Import_Stub, pstate, imp->incs()[i]));
      }
    }

    else if (lex < kwd_extend >(true)) {
      Lookahead lookahead = lookahead_for_include(position);
      if (!lookahead.found) css_error("Invalid CSS", " after ", ": expected selector, was ");
      Selector_List_Obj target;
      if (!lookahead.has_interpolants) {
        target = parse_selector_list(true);
      }
      else {
        target = SASS_MEMORY_NEW(Selector_List, pstate);
        target->schema(parse_selector_schema(lookahead.found, true));
      }

      block->append(SASS_MEMORY_NEW(Extension, pstate, target));
    }

    // selector may contain interpolations which need delayed evaluation
    else if (
      !(lookahead_result = lookahead_for_selector(position)).error &&
      !lookahead_result.is_custom_property
    )
    {
      block->append(parse_ruleset(lookahead_result));
    }

    // parse multiple specific keyword directives
    else if (lex < kwd_media >(true)) { block->append(parse_media_block()); }
    else if (lex < kwd_at_root >(true)) { block->append(parse_at_root_block()); }
    else if (lex < kwd_include_directive >(true)) { block->append(parse_include_directive()); }
    else if (lex < kwd_content_directive >(true)) { block->append(parse_content_directive()); }
    else if (lex < kwd_supports_directive >(true)) { block->append(parse_supports_directive()); }
    else if (lex < kwd_mixin >(true)) { block->append(parse_definition(Definition::MIXIN)); }
    else if (lex < kwd_function >(true)) { block->append(parse_definition(Definition::FUNCTION)); }

    // ignore the @charset directive for now
    else if (lex< kwd_charset_directive >(true)) { parse_charset_directive(); }

    // generic at keyword (keep last)
    else if (lex< re_special_directive >(true)) { block->append(parse_special_directive()); }
    else if (lex< re_prefixed_directive >(true)) { block->append(parse_prefixed_directive()); }
    else if (lex< at_keyword >(true)) { block->append(parse_directive()); }

    else if (is_root && stack.back() != Scope::AtRoot /* && block->is_root() */) {
      lex< css_whitespace >();
      if (position >= end) return true;
      css_error("Invalid CSS", " after ", ": expected 1 selector or at-rule, was ");
    }
    // parse a declaration
    else
    {
      // ToDo: how does it handle parse errors?
      // maybe we are expected to parse something?
      Declaration_Obj decl = parse_declaration();
      decl->tabs(indentation);
      block->append(decl);
      // maybe we have a "sub-block"
      if (peek< exactly<'{'> >()) {
        if (decl->is_indented()) ++ indentation;
        // parse a propset that rides on the declaration's property
        stack.push_back(Scope::Properties);
        decl->block(parse_block());
        stack.pop_back();
        if (decl->is_indented()) -- indentation;
      }
    }
    // something matched
    return true;
  }
  // EO parse_block_nodes

  // parse imports inside the
  Import_Obj Parser::parse_import()
  {
    Import_Obj imp = SASS_MEMORY_NEW(Import, pstate);
    std::vector<std::pair<std::string,Function_Call_Obj>> to_import;
    bool first = true;
    do {
      while (lex< block_comment >());
      if (lex< quoted_string >()) {
        to_import.push_back(std::pair<std::string,Function_Call_Obj>(std::string(lexed), {}));
      }
      else if (lex< uri_prefix >()) {
        Arguments_Obj args = SASS_MEMORY_NEW(Arguments, pstate);
        Function_Call_Obj result = SASS_MEMORY_NEW(Function_Call, pstate, std::string("url"), args);

        if (lex< quoted_string >()) {
          Expression_Obj quoted_url = parse_string();
          args->append(SASS_MEMORY_NEW(Argument, quoted_url->pstate(), quoted_url));
        }
        else if (String_Obj string_url = parse_url_function_argument()) {
          args->append(SASS_MEMORY_NEW(Argument, string_url->pstate(), string_url));
        }
        else if (peek < skip_over_scopes < exactly < '(' >, exactly < ')' > > >(position)) {
          Expression_Obj braced_url = parse_list(); // parse_interpolated_chunk(lexed);
          args->append(SASS_MEMORY_NEW(Argument, braced_url->pstate(), braced_url));
        }
        else {
          error("malformed URL");
        }
        if (!lex< exactly<')'> >()) error("URI is missing ')'");
        to_import.push_back(std::pair<std::string, Function_Call_Obj>("", result));
      }
      else {
        if (first) error("@import directive requires a url or quoted path");
        else error("expecting another url or quoted path in @import list");
      }
      first = false;
    } while (lex_css< exactly<','> >());

    if (!peek_css< alternatives< exactly<';'>, exactly<'}'>, end_of_file > >()) {
      List_Obj import_queries = parse_media_queries();
      imp->import_queries(import_queries);
    }

    for(auto location : to_import) {
      if (location.second) {
        imp->urls().push_back(location.second);
      }
      // check if custom importers want to take over the handling
      else if (!ctx.call_importers(unquote(location.first), path, pstate, imp)) {
        // nobody wants it, so we do our import
        ctx.import_url(imp, location.first, path);
      }
    }

    return imp;
  }

  Definition_Obj Parser::parse_definition(Definition::Type which_type)
  {
    std::string which_str(lexed);
    if (!lex< identifier >()) error("invalid name in " + which_str + " definition");
    std::string name(Util::normalize_underscores(lexed));
    if (which_type == Definition::FUNCTION && (name == "and" || name == "or" || name == "not"))
    { error("Invalid function name \"" + name + "\"."); }
    ParserState source_position_of_def = pstate;
    Parameters_Obj params = parse_parameters();
    if (which_type == Definition::MIXIN) stack.push_back(Scope::Mixin);
    else stack.push_back(Scope::Function);
    Block_Obj body = parse_block();
    stack.pop_back();
    return SASS_MEMORY_NEW(Definition, source_position_of_def, name, params, body, which_type);
  }

  Parameters_Obj Parser::parse_parameters()
  {
    Parameters_Obj params = SASS_MEMORY_NEW(Parameters, pstate);
    if (lex_css< exactly<'('> >()) {
      // if there's anything there at all
      if (!peek_css< exactly<')'> >()) {
        do {
          if (peek< exactly<')'> >()) break;
          params->append(parse_parameter());
        } while (lex_css< exactly<','> >());
      }
      if (!lex_css< exactly<')'> >()) {
        css_error("Invalid CSS", " after ", ": expected \")\", was ");
      }
    }
    return params;
  }

  Parameter_Obj Parser::parse_parameter()
  {
    if (peek< alternatives< exactly<','>, exactly< '{' >, exactly<';'> > >()) {
      css_error("Invalid CSS", " after ", ": expected variable (e.g. $foo), was ");
    }
    while (lex< alternatives < spaces, block_comment > >());
    lex < variable >();
    std::string name(Util::normalize_underscores(lexed));
    ParserState pos = pstate;
    Expression_Obj val;
    bool is_rest = false;
    while (lex< alternatives < spaces, block_comment > >());
    if (lex< exactly<':'> >()) { // there's a default value
      while (lex< block_comment >());
      val = parse_space_list();
    }
    else if (lex< exactly< ellipsis > >()) {
      is_rest = true;
    }
    return SASS_MEMORY_NEW(Parameter, pos, name, val, is_rest);
  }

  Arguments_Obj Parser::parse_arguments()
  {
    Arguments_Obj args = SASS_MEMORY_NEW(Arguments, pstate);
    if (lex_css< exactly<'('> >()) {
      // if there's anything there at all
      if (!peek_css< exactly<')'> >()) {
        do {
          if (peek< exactly<')'> >()) break;
          args->append(parse_argument());
        } while (lex_css< exactly<','> >());
      }
      if (!lex_css< exactly<')'> >()) {
        css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
      }
    }
    return args;
  }

  Argument_Obj Parser::parse_argument()
  {
    if (peek< alternatives< exactly<','>, exactly< '{' >, exactly<';'> > >()) {
      css_error("Invalid CSS", " after ", ": expected \")\", was ");
    }
    if (peek_css< sequence < exactly< hash_lbrace >, exactly< rbrace > > >()) {
      position += 2;
      css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
    }

    Argument_Obj arg;
    if (peek_css< sequence < variable, optional_css_comments, exactly<':'> > >()) {
      lex_css< variable >();
      std::string name(Util::normalize_underscores(lexed));
      ParserState p = pstate;
      lex_css< exactly<':'> >();
      Expression_Obj val = parse_space_list();
      arg = SASS_MEMORY_NEW(Argument, p, val, name);
    }
    else {
      bool is_arglist = false;
      bool is_keyword = false;
      Expression_Obj val = parse_space_list();
      List_Ptr l = Cast<List>(val);
      if (lex_css< exactly< ellipsis > >()) {
        if (val->concrete_type() == Expression::MAP || (
           (l != NULL && l->separator() == SASS_HASH)
        )) is_keyword = true;
        else is_arglist = true;
      }
      arg = SASS_MEMORY_NEW(Argument, pstate, val, "", is_arglist, is_keyword);
    }
    return arg;
  }

  Assignment_Obj Parser::parse_assignment()
  {
    std::string name(Util::normalize_underscores(lexed));
    ParserState var_source_position = pstate;
    if (!lex< exactly<':'> >()) error("expected ':' after " + name + " in assignment statement");
    if (peek_css< alternatives < exactly<';'>, end_of_file > >()) {
      css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
    }
    Expression_Obj val;
    Lookahead lookahead = lookahead_for_value(position);
    if (lookahead.has_interpolants && lookahead.found) {
      val = parse_value_schema(lookahead.found);
    } else {
      val = parse_list();
    }
    bool is_default = false;
    bool is_global = false;
    while (peek< alternatives < default_flag, global_flag > >()) {
      if (lex< default_flag >()) is_default = true;
      else if (lex< global_flag >()) is_global = true;
    }
    return SASS_MEMORY_NEW(Assignment, var_source_position, name, val, is_default, is_global);
  }

  // a ruleset connects a selector and a block
  Ruleset_Obj Parser::parse_ruleset(Lookahead lookahead)
  {
    NESTING_GUARD(nestings);
    // inherit is_root from parent block
    Block_Obj parent = block_stack.back();
    bool is_root = parent && parent->is_root();
    // make sure to move up the the last position
    lex < optional_css_whitespace >(false, true);
    // create the connector object (add parts later)
    Ruleset_Obj ruleset = SASS_MEMORY_NEW(Ruleset, pstate);
    // parse selector static or as schema to be evaluated later
    if (lookahead.parsable) ruleset->selector(parse_selector_list(false));
    else {
      Selector_List_Obj list = SASS_MEMORY_NEW(Selector_List, pstate);
      list->schema(parse_selector_schema(lookahead.position, false));
      ruleset->selector(list);
    }
    // then parse the inner block
    stack.push_back(Scope::Rules);
    ruleset->block(parse_block());
    stack.pop_back();
    // update for end position
    ruleset->update_pstate(pstate);
    ruleset->block()->update_pstate(pstate);
    // need this info for sanity checks
    ruleset->is_root(is_root);
    // return AST Node
    return ruleset;
  }

  // parse a selector schema that will be evaluated in the eval stage
  // uses a string schema internally to do the actual schema handling
  // in the eval stage we will be re-parse it into an actual selector
  Selector_Schema_Obj Parser::parse_selector_schema(const char* end_of_selector, bool chroot)
  {
    NESTING_GUARD(nestings);
    // move up to the start
    lex< optional_spaces >();
    const char* i = position;
    // selector schema re-uses string schema implementation
    String_Schema_Ptr schema = SASS_MEMORY_NEW(String_Schema, pstate);
    // the selector schema is pretty much just a wrapper for the string schema
    Selector_Schema_Obj selector_schema = SASS_MEMORY_NEW(Selector_Schema, pstate, schema);
    selector_schema->connect_parent(chroot == false);
    selector_schema->media_block(last_media_block);

    // process until end
    while (i < end_of_selector) {
      // try to parse mutliple interpolants
      if (const char* p = find_first_in_interval< exactly<hash_lbrace>, block_comment >(i, end_of_selector)) {
        // accumulate the preceding segment if the position has advanced
        if (i < p) {
          std::string parsed(i, p);
          String_Constant_Obj str = SASS_MEMORY_NEW(String_Constant, pstate, parsed);
          pstate += Offset(parsed);
          str->update_pstate(pstate);
          schema->append(str);
        }

        // skip over all nested inner interpolations up to our own delimiter
        const char* j = skip_over_scopes< exactly<hash_lbrace>, exactly<rbrace> >(p + 2, end_of_selector);
        // check if the interpolation never ends of only contains white-space (error out)
        if (!j || peek < sequence < optional_spaces, exactly<rbrace> > >(p+2)) {
          position = p+2;
          css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
        }
        // pass inner expression to the parser to resolve nested interpolations
        pstate.add(p, p+2);
        Expression_Obj interpolant = Parser::from_c_str(p+2, j, ctx, traces, pstate).parse_list();
        // set status on the list expression
        interpolant->is_interpolant(true);
        // schema->has_interpolants(true);
        // add to the string schema
        schema->append(interpolant);
        // advance parser state
        pstate.add(p+2, j);
        // advance position
        i = j;
      }
      // no more interpolants have been found
      // add the last segment if there is one
      else {
        // make sure to add the last bits of the string up to the end (if any)
        if (i < end_of_selector) {
          std::string parsed(i, end_of_selector);
          String_Constant_Obj str = SASS_MEMORY_NEW(String_Constant, pstate, parsed);
          pstate += Offset(parsed);
          str->update_pstate(pstate);
          i = end_of_selector;
          schema->append(str);
        }
        // exit loop
      }
    }
    // EO until eos

    // update position
    position = i;

    // update for end position
    selector_schema->update_pstate(pstate);
    schema->update_pstate(pstate);

    after_token = before_token = pstate;

    // return parsed result
    return selector_schema.detach();
  }
  // EO parse_selector_schema

  void Parser::parse_charset_directive()
  {
    lex <
      sequence <
        quoted_string,
        optional_spaces,
        exactly <';'>
      >
    >();
  }

  // called after parsing `kwd_include_directive`
  Mixin_Call_Obj Parser::parse_include_directive()
  {
    // lex identifier into `lexed` var
    lex_identifier(); // may error out
    // normalize underscores to hyphens
    std::string name(Util::normalize_underscores(lexed));
    // create the initial mixin call object
    Mixin_Call_Obj call = SASS_MEMORY_NEW(Mixin_Call, pstate, name, {}, {}, {});
    // parse mandatory arguments
    call->arguments(parse_arguments());
    // parse using and optional block parameters
    bool has_parameters = lex< kwd_using >() != nullptr;

    if (has_parameters) {
      if (!peek< exactly<'('> >()) css_error("Invalid CSS", " after ", ": expected \"(\", was ");
    } else {
      if (peek< exactly<'('> >()) css_error("Invalid CSS", " after ", ": expected \";\", was ");
    }

    if (has_parameters) call->block_parameters(parse_parameters());

    // parse optional block
    if (peek < exactly <'{'> >()) {
      call->block(parse_block());
    }
    else if (has_parameters)  {
      css_error("Invalid CSS", " after ", ": expected \"{\", was ");
    }
    // return ast node
    return call.detach();
  }
  // EO parse_include_directive

  // parse a list of complex selectors
  // this is the main entry point for most
  Selector_List_Obj Parser::parse_selector_list(bool chroot)
  {
    bool reloop;
    bool had_linefeed = false;
    NESTING_GUARD(nestings);
    Complex_Selector_Obj sel;
    Selector_List_Obj group = SASS_MEMORY_NEW(Selector_List, pstate);
    group->media_block(last_media_block);

    if (peek_css< alternatives < end_of_file, exactly <'{'>, exactly <','> > >()) {
      css_error("Invalid CSS", " after ", ": expected selector, was ");
    }

    do {
      reloop = false;

      had_linefeed = had_linefeed || peek_newline();

      if (peek_css< alternatives < class_char < selector_list_delims > > >())
        break; // in case there are superfluous commas at the end

      // now parse the complex selector
      sel = parse_complex_selector(chroot);

      if (!sel) return group.detach();

      sel->has_line_feed(had_linefeed);

      had_linefeed = false;

      while (peek_css< exactly<','> >())
      {
        lex< css_comments >(false);
        // consume everything up and including the comma separator
        reloop = lex< exactly<','> >() != 0;
        // remember line break (also between some commas)
        had_linefeed = had_linefeed || peek_newline();
        // remember line break (also between some commas)
      }
      group->append(sel);
    }
    while (reloop);
    while (lex_css< kwd_optional >()) {
      group->is_optional(true);
    }
    // update for end position
    group->update_pstate(pstate);
    if (sel) sel->mutable_last()->has_line_break(false);
    return group.detach();
  }
  // EO parse_selector_list

  // a complex selector combines a compound selector with another
  // complex selector, with one of four combinator operations.
  // the compound selector (head) is optional, since the combinator
  // can come first in the whole selector sequence (like `> DIV').
  Complex_Selector_Obj Parser::parse_complex_selector(bool chroot)
  {

    NESTING_GUARD(nestings);
    String_Obj reference;
    lex < block_comment >();
    advanceToNextToken();
    Complex_Selector_Obj sel = SASS_MEMORY_NEW(Complex_Selector, pstate);

    if (peek < end_of_file >()) return {};

    // parse the left hand side
    Compound_Selector_Obj lhs;
    // special case if it starts with combinator ([+~>])
    if (!peek_css< class_char < selector_combinator_ops > >()) {
      // parse the left hand side
      lhs = parse_compound_selector();
    }


    // parse combinator between lhs and rhs
    Complex_Selector::Combinator combinator = Complex_Selector::ANCESTOR_OF;
    if      (lex< exactly<'+'> >()) combinator = Complex_Selector::ADJACENT_TO;
    else if (lex< exactly<'~'> >()) combinator = Complex_Selector::PRECEDES;
    else if (lex< exactly<'>'> >()) combinator = Complex_Selector::PARENT_OF;
    else if (lex< sequence < exactly<'/'>, negate < exactly < '*' > > > >()) {
      // comments are allowed, but not spaces?
      combinator = Complex_Selector::REFERENCE;
      if (!lex < re_reference_combinator >()) return {};
      reference = SASS_MEMORY_NEW(String_Constant, pstate, lexed);
      if (!lex < exactly < '/' > >()) return {}; // ToDo: error msg?
    }

    if (!lhs && combinator == Complex_Selector::ANCESTOR_OF) return {};

    // lex < block_comment >();
    sel->head(lhs);
    sel->combinator(combinator);
    sel->media_block(last_media_block);

    if (combinator == Complex_Selector::REFERENCE) sel->reference(reference);
    // has linfeed after combinator?
    sel->has_line_break(peek_newline());
    // sel->has_line_feed(has_line_feed);

    // check if we got the abort condition (ToDo: optimize)
    if (!peek_css< class_char < complex_selector_delims > >()) {
      // parse next selector in sequence
      sel->tail(parse_complex_selector(true));
    }

    // add a parent selector if we are not in a root
    // also skip adding parent ref if we only have refs
    if (!sel->has_parent_ref() && !chroot) {
      // create the objects to wrap parent selector reference
      Compound_Selector_Obj head = SASS_MEMORY_NEW(Compound_Selector, pstate);
      Parent_Selector_Ptr parent = SASS_MEMORY_NEW(Parent_Selector, pstate, false);
      parent->media_block(last_media_block);
      head->media_block(last_media_block);
      // add simple selector
      head->append(parent);
      // selector may not have any head yet
      if (!sel->head()) { sel->head(head); }
      // otherwise we need to create a new complex selector and set the old one as its tail
      else {
        sel = SASS_MEMORY_NEW(Complex_Selector, pstate, Complex_Selector::ANCESTOR_OF, head, sel);
        sel->media_block(last_media_block);
      }
      // peek for linefeed and remember result on head
      // if (peek_newline()) head->has_line_break(true);
    }

    sel->update_pstate(pstate);
    // complex selector
    return sel;
  }
  // EO parse_complex_selector

  // parse one compound selector, which is basically
  // a list of simple selectors (directly adjacent)
  // lex them exactly (without skipping white-space)
  Compound_Selector_Obj Parser::parse_compound_selector()
  {
    // init an empty compound selector wrapper
    Compound_Selector_Obj seq = SASS_MEMORY_NEW(Compound_Selector, pstate);
    seq->media_block(last_media_block);

    // skip initial white-space
    lex< css_whitespace >();

    // parse list
    while (true)
    {
      // remove all block comments (don't skip white-space)
      lex< delimited_by< slash_star, star_slash, false > >(false);
      // parse functional
      if (match < re_pseudo_selector >())
      {
        seq->append(parse_simple_selector());
      }
      // parse parent selector
      else if (lex< exactly<'&'> >(false))
      {
        if (!allow_parent) error("Parent selectors aren't allowed here.");
        // this produces a linefeed!?
        seq->has_parent_reference(true);
        seq->append(SASS_MEMORY_NEW(Parent_Selector, pstate));
        // parent selector only allowed at start
        // upcoming Sass may allow also trailing
        if (seq->length() > 1) {
          ParserState state(pstate);
          Simple_Selector_Obj cur = (*seq)[seq->length()-1];
          Simple_Selector_Obj prev = (*seq)[seq->length()-2];
          std::string sel(prev->to_string({ NESTED, 5 }));
          std::string found(cur->to_string({ NESTED, 5 }));
          if (lex < identifier >()) { found += std::string(lexed); }
          error("Invalid CSS after \"" + sel + "\": expected \"{\", was \"" + found + "\"\n\n"
            "\"" + found + "\" may only be used at the beginning of a compound selector.", state);
        }
      }
      // parse type selector
      else if (lex< re_type_selector >(false))
      {
        seq->append(SASS_MEMORY_NEW(Type_Selector, pstate, lexed));
      }
      // peek for abort conditions
      else if (peek< spaces >()) break;
      else if (peek< end_of_file >()) { break; }
      else if (peek_css < class_char < selector_combinator_ops > >()) break;
      else if (peek_css < class_char < complex_selector_delims > >()) break;
      // otherwise parse another simple selector
      else {
        Simple_Selector_Obj sel = parse_simple_selector();
        if (!sel) return {};
        seq->append(sel);
      }
    }

    if (seq && !peek_css<alternatives<end_of_file,exactly<'{'>>>()) {
      seq->has_line_break(peek_newline());
    }

    // EO while true
    return seq;

  }
  // EO parse_compound_selector

  Simple_Selector_Obj Parser::parse_simple_selector()
  {
    lex < css_comments >(false);
    if (lex< class_name >()) {
      return SASS_MEMORY_NEW(Class_Selector, pstate, lexed);
    }
    else if (lex< id_name >()) {
      return SASS_MEMORY_NEW(Id_Selector, pstate, lexed);
    }
    else if (lex< alternatives < variable, number, static_reference_combinator > >()) {
      return SASS_MEMORY_NEW(Type_Selector, pstate, lexed);
    }
    else if (peek< pseudo_not >()) {
      return parse_negated_selector();
    }
    else if (peek< re_pseudo_selector >()) {
      return parse_pseudo_selector();
    }
    else if (peek< exactly<':'> >()) {
      return parse_pseudo_selector();
    }
    else if (lex < exactly<'['> >()) {
      return parse_attribute_selector();
    }
    else if (lex< placeholder >()) {
      Placeholder_Selector_Ptr sel = SASS_MEMORY_NEW(Placeholder_Selector, pstate, lexed);
      sel->media_block(last_media_block);
      return sel;
    }
    else {
      css_error("Invalid CSS", " after ", ": expected selector, was ");
    }
    // failed
    return {};
  }

  Wrapped_Selector_Obj Parser::parse_negated_selector()
  {
    lex< pseudo_not >();
    std::string name(lexed);
    ParserState nsource_position = pstate;
    Selector_List_Obj negated = parse_selector_list(true);
    if (!lex< exactly<')'> >()) {
      error("negated selector is missing ')'");
    }
    name.erase(name.size() - 1);
    return SASS_MEMORY_NEW(Wrapped_Selector, nsource_position, name, negated);
  }

  // a pseudo selector often starts with one or two colons
  // it can contain more selectors inside parentheses
  Simple_Selector_Obj Parser::parse_pseudo_selector() {
    if (lex< sequence<
          pseudo_prefix,
          // we keep the space within the name, strange enough
          // ToDo: refactor output to schedule the space for it
          // or do we really want to keep the real white-space?
          sequence< identifier, optional < block_comment >, exactly<'('> >
        > >())
    {

      std::string name(lexed);
      name.erase(name.size() - 1);
      ParserState p = pstate;

      // specially parse static stuff
      // ToDo: really everything static?
      if (peek_css <
            sequence <
              alternatives <
                static_value,
                binomial
              >,
              optional_css_whitespace,
              exactly<')'>
            >
          >()
      ) {
        lex_css< alternatives < static_value, binomial > >();
        String_Constant_Obj expr = SASS_MEMORY_NEW(String_Constant, pstate, lexed);
        if (lex_css< exactly<')'> >()) {
          expr->can_compress_whitespace(true);
          return SASS_MEMORY_NEW(Pseudo_Selector, p, name, expr);
        }
      }
      else if (Selector_List_Obj wrapped = parse_selector_list(true)) {
        if (wrapped && lex_css< exactly<')'> >()) {
          return SASS_MEMORY_NEW(Wrapped_Selector, p, name, wrapped);
        }
      }

    }
    // EO if pseudo selector

    else if (lex < sequence< optional < pseudo_prefix >, identifier > >()) {
      return SASS_MEMORY_NEW(Pseudo_Selector, pstate, lexed);
    }
    else if(lex < pseudo_prefix >()) {
      css_error("Invalid CSS", " after ", ": expected pseudoclass or pseudoelement, was ");
    }

    css_error("Invalid CSS", " after ", ": expected \")\", was ");

    // unreachable statement
    return {};
  }

  const char* Parser::re_attr_sensitive_close(const char* src)
  {
    return alternatives < exactly<']'>, exactly<'/'> >(src);
  }

  const char* Parser::re_attr_insensitive_close(const char* src)
  {
    return sequence < insensitive<'i'>, re_attr_sensitive_close >(src);
  }

  Attribute_Selector_Obj Parser::parse_attribute_selector()
  {
    ParserState p = pstate;
    if (!lex_css< attribute_name >()) error("invalid attribute name in attribute selector");
    std::string name(lexed);
    if (lex_css< re_attr_sensitive_close >()) {
      return SASS_MEMORY_NEW(Attribute_Selector, p, name, "", {}, {});
    }
    else if (lex_css< re_attr_insensitive_close >()) {
      char modifier = lexed.begin[0];
      return SASS_MEMORY_NEW(Attribute_Selector, p, name, "", {}, modifier);
    }
    if (!lex_css< alternatives< exact_match, class_match, dash_match,
                                prefix_match, suffix_match, substring_match > >()) {
      error("invalid operator in attribute selector for " + name);
    }
    std::string matcher(lexed);

    String_Obj value;
    if (lex_css< identifier >()) {
      value = SASS_MEMORY_NEW(String_Constant, p, lexed);
    }
    else if (lex_css< quoted_string >()) {
      value = parse_interpolated_chunk(lexed, true); // needed!
    }
    else {
      error("expected a string constant or identifier in attribute selector for " + name);
    }

    if (lex_css< re_attr_sensitive_close >()) {
      return SASS_MEMORY_NEW(Attribute_Selector, p, name, matcher, value, 0);
    }
    else if (lex_css< re_attr_insensitive_close >()) {
      char modifier = lexed.begin[0];
      return SASS_MEMORY_NEW(Attribute_Selector, p, name, matcher, value, modifier);
    }
    error("unterminated attribute selector for " + name);
    return {}; // to satisfy compilers (error must not return)
  }

  /* parse block comment and add to block */
  void Parser::parse_block_comments()
  {
    Block_Obj block = block_stack.back();

    while (lex< block_comment >()) {
      bool is_important = lexed.begin[2] == '!';
      // flag on second param is to skip loosely over comments
      String_Obj contents = parse_interpolated_chunk(lexed, true, false);
      block->append(SASS_MEMORY_NEW(Comment, pstate, contents, is_important));
    }
  }

  Declaration_Obj Parser::parse_declaration() {
    String_Obj prop;
    bool is_custom_property = false;
    if (lex< sequence< optional< exactly<'*'> >, identifier_schema > >()) {
      const std::string property(lexed);
      is_custom_property = property.compare(0, 2, "--") == 0;
      prop = parse_identifier_schema();
    }
    else if (lex< sequence< optional< exactly<'*'> >, identifier, zero_plus< block_comment > > >()) {
      const std::string property(lexed);
      is_custom_property = property.compare(0, 2, "--") == 0;
      prop = SASS_MEMORY_NEW(String_Constant, pstate, lexed);
    }
    else {
      css_error("Invalid CSS", " after ", ": expected \"}\", was ");
    }
    bool is_indented = true;
    const std::string property(lexed);
    if (!lex_css< one_plus< exactly<':'> > >()) error("property \"" + escape_string(property)  + "\" must be followed by a ':'");
    if (!is_custom_property && match< sequence< optional_css_comments, exactly<';'> > >()) error("style declaration must contain a value");
    if (match< sequence< optional_css_comments, exactly<'{'> > >()) is_indented = false; // don't indent if value is empty
    if (is_custom_property) {
      return SASS_MEMORY_NEW(Declaration, prop->pstate(), prop, parse_css_variable_value(), false, true);
    }
    lex < css_comments >(false);
    if (peek_css< static_value >()) {
      return SASS_MEMORY_NEW(Declaration, prop->pstate(), prop, parse_static_value()/*, lex<kwd_important>()*/);
    }
    else {
      Expression_Obj value;
      Lookahead lookahead = lookahead_for_value(position);
      if (lookahead.found) {
        if (lookahead.has_interpolants) {
          value = parse_value_schema(lookahead.found);
        } else {
          value = parse_list(DELAYED);
        }
      }
      else {
        value = parse_list(DELAYED);
        if (List_Ptr list = Cast<List>(value)) {
          if (!list->is_bracketed() && list->length() == 0 && !peek< exactly <'{'> >()) {
            css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
          }
        }
      }
      lex < css_comments >(false);
      Declaration_Obj decl = SASS_MEMORY_NEW(Declaration, prop->pstate(), prop, value/*, lex<kwd_important>()*/);
      decl->is_indented(is_indented);
      decl->update_pstate(pstate);
      return decl;
    }
  }

  // parse +/- and return false if negative
  // this is never hit via spec tests
  bool Parser::parse_number_prefix()
  {
    bool positive = true;
    while(true) {
      if (lex < block_comment >()) continue;
      if (lex < number_prefix >()) continue;
      if (lex < exactly < '-' > >()) {
        positive = !positive;
        continue;
      }
      break;
    }
    return positive;
  }

  Expression_Obj Parser::parse_map()
  {
    NESTING_GUARD(nestings);
    Expression_Obj key = parse_list();
    List_Obj map = SASS_MEMORY_NEW(List, pstate, 0, SASS_HASH);

    // it's not a map so return the lexed value as a list value
    if (!lex_css< exactly<':'> >())
    { return key; }

    List_Obj l = Cast<List>(key);
    if (l && l->separator() == SASS_COMMA) {
      css_error("Invalid CSS", " after ", ": expected \")\", was ");
    }

    Expression_Obj value = parse_space_list();

    map->append(key);
    map->append(value);

    while (lex_css< exactly<','> >())
    {
      // allow trailing commas - #495
      if (peek_css< exactly<')'> >(position))
      { break; }

      key = parse_space_list();

      if (!(lex< exactly<':'> >()))
      { css_error("Invalid CSS", " after ", ": expected \":\", was "); }

      value = parse_space_list();

      map->append(key);
      map->append(value);
    }

    ParserState ps = map->pstate();
    ps.offset = pstate - ps + pstate.offset;
    map->pstate(ps);

    return map;
  }

  Expression_Obj Parser::parse_bracket_list()
  {
    NESTING_GUARD(nestings);
    // check if we have an empty list
    // return the empty list as such
    if (peek_css< list_terminator >(position))
    {
      // return an empty list (nothing to delay)
      return SASS_MEMORY_NEW(List, pstate, 0, SASS_SPACE, false, true);
    }

    bool has_paren = peek_css< exactly<'('> >() != NULL;

    // now try to parse a space list
    Expression_Obj list = parse_space_list();
    // if it's a singleton, return it (don't wrap it)
    if (!peek_css< exactly<','> >(position)) {
      List_Obj l = Cast<List>(list);
      if (!l || l->is_bracketed() || has_paren) {
        List_Obj bracketed_list = SASS_MEMORY_NEW(List, pstate, 1, SASS_SPACE, false, true);
        bracketed_list->append(list);
        return bracketed_list;
      }
      l->is_bracketed(true);
      return l;
    }

    // if we got so far, we actually do have a comma list
    List_Obj bracketed_list = SASS_MEMORY_NEW(List, pstate, 2, SASS_COMMA, false, true);
    // wrap the first expression
    bracketed_list->append(list);

    while (lex_css< exactly<','> >())
    {
      // check for abort condition
      if (peek_css< list_terminator >(position)
      ) { break; }
      // otherwise add another expression
      bracketed_list->append(parse_space_list());
    }
    // return the list
    return bracketed_list;
  }

  // parse list returns either a space separated list,
  // a comma separated list or any bare expression found.
  // so to speak: we unwrap items from lists if possible here!
  Expression_Obj Parser::parse_list(bool delayed)
  {
    NESTING_GUARD(nestings);
    return parse_comma_list(delayed);
  }

  // will return singletons unwrapped
  Expression_Obj Parser::parse_comma_list(bool delayed)
  {
    NESTING_GUARD(nestings);
    // check if we have an empty list
    // return the empty list as such
    if (peek_css< list_terminator >(position))
    {
      // return an empty list (nothing to delay)
      return SASS_MEMORY_NEW(List, pstate, 0);
    }

    // now try to parse a space list
    Expression_Obj list = parse_space_list();
    // if it's a singleton, return it (don't wrap it)
    if (!peek_css< exactly<','> >(position)) {
      // set_delay doesn't apply to list children
      // so this will only undelay single values
      if (!delayed) list->set_delayed(false);
      return list;
    }

    // if we got so far, we actually do have a comma list
    List_Obj comma_list = SASS_MEMORY_NEW(List, pstate, 2, SASS_COMMA);
    // wrap the first expression
    comma_list->append(list);

    while (lex_css< exactly<','> >())
    {
      // check for abort condition
      if (peek_css< list_terminator >(position)
      ) { break; }
      // otherwise add another expression
      comma_list->append(parse_space_list());
    }
    // return the list
    return comma_list;
  }
  // EO parse_comma_list

  // will return singletons unwrapped
  Expression_Obj Parser::parse_space_list()
  {
    NESTING_GUARD(nestings);
    Expression_Obj disj1 = parse_disjunction();
    // if it's a singleton, return it (don't wrap it)
    if (peek_css< space_list_terminator >(position)
    ) {
      return disj1; }

    List_Obj space_list = SASS_MEMORY_NEW(List, pstate, 2, SASS_SPACE);
    space_list->append(disj1);

    while (
      !(peek_css< space_list_terminator >(position)) &&
      peek_css< optional_css_whitespace >() != end
    ) {
      // the space is parsed implicitly?
      space_list->append(parse_disjunction());
    }
    // return the list
    return space_list;
  }
  // EO parse_space_list

  // parse logical OR operation
  Expression_Obj Parser::parse_disjunction()
  {
    NESTING_GUARD(nestings);
    advanceToNextToken();
    ParserState state(pstate);
    // parse the left hand side conjunction
    Expression_Obj conj = parse_conjunction();
    // parse multiple right hand sides
    std::vector<Expression_Obj> operands;
    while (lex_css< kwd_or >())
      operands.push_back(parse_conjunction());
    // if it's a singleton, return it directly
    if (operands.size() == 0) return conj;
    // fold all operands into one binary expression
    Expression_Obj ex = fold_operands(conj, operands, { Sass_OP::OR });
    state.offset = pstate - state + pstate.offset;
    ex->pstate(state);
    return ex;
  }
  // EO parse_disjunction

  // parse logical AND operation
  Expression_Obj Parser::parse_conjunction()
  {
    NESTING_GUARD(nestings);
    advanceToNextToken();
    ParserState state(pstate);
    // parse the left hand side relation
    Expression_Obj rel = parse_relation();
    // parse multiple right hand sides
    std::vector<Expression_Obj> operands;
    while (lex_css< kwd_and >()) {
      operands.push_back(parse_relation());
    }
    // if it's a singleton, return it directly
    if (operands.size() == 0) return rel;
    // fold all operands into one binary expression
    Expression_Obj ex = fold_operands(rel, operands, { Sass_OP::AND });
    state.offset = pstate - state + pstate.offset;
    ex->pstate(state);
    return ex;
  }
  // EO parse_conjunction

  // parse comparison operations
  Expression_Obj Parser::parse_relation()
  {
    NESTING_GUARD(nestings);
    advanceToNextToken();
    ParserState state(pstate);
    // parse the left hand side expression
    Expression_Obj lhs = parse_expression();
    std::vector<Expression_Obj> operands;
    std::vector<Operand> operators;
    // if it's a singleton, return it (don't wrap it)
    while (peek< alternatives <
            kwd_eq,
            kwd_neq,
            kwd_gte,
            kwd_gt,
            kwd_lte,
            kwd_lt
          > >(position))
    {
      // is directly adjancent to expression?
      bool left_ws = peek < css_comments >() != NULL;
      // parse the operator
      enum Sass_OP op
      = lex<kwd_eq>()  ? Sass_OP::EQ
      : lex<kwd_neq>() ? Sass_OP::NEQ
      : lex<kwd_gte>() ? Sass_OP::GTE
      : lex<kwd_lte>() ? Sass_OP::LTE
      : lex<kwd_gt>()  ? Sass_OP::GT
      : lex<kwd_lt>()  ? Sass_OP::LT
      // we checked the possibilities on top of fn
      :                  Sass_OP::EQ;
      // is directly adjacent to expression?
      bool right_ws = peek < css_comments >() != NULL;
      operators.push_back({ op, left_ws, right_ws });
      operands.push_back(parse_expression());
    }
    // we are called recursively for list, so we first
    // fold inner binary expression which has delayed
    // correctly set to zero. After folding we also unwrap
    // single nested items. So we cannot set delay on the
    // returned result here, as we have lost nestings ...
    Expression_Obj ex = fold_operands(lhs, operands, operators);
    state.offset = pstate - state + pstate.offset;
    ex->pstate(state);
    return ex;
  }
  // parse_relation

  // parse expression valid for operations
  // called from parse_relation
  // called from parse_for_directive
  // called from parse_media_expression
  // parse addition and subtraction operations
  Expression_Obj Parser::parse_expression()
  {
    NESTING_GUARD(nestings);
    advanceToNextToken();
    ParserState state(pstate);
    // parses multiple add and subtract operations
    // NOTE: make sure that identifiers starting with
    // NOTE: dashes do NOT count as subtract operation
    Expression_Obj lhs = parse_operators();
    // if it's a singleton, return it (don't wrap it)
    if (!(peek_css< exactly<'+'> >(position) ||
          // condition is a bit misterious, but some combinations should not be counted as operations
          (peek< no_spaces >(position) && peek< sequence< negate< unsigned_number >, exactly<'-'>, negate< space > > >(position)) ||
          (peek< sequence< negate< unsigned_number >, exactly<'-'>, negate< unsigned_number > > >(position))) ||
          peek< sequence < zero_plus < exactly <'-' > >, identifier > >(position))
    { return lhs; }

    std::vector<Expression_Obj> operands;
    std::vector<Operand> operators;
    bool left_ws = peek < css_comments >() != NULL;
    while (
      lex_css< exactly<'+'> >() ||

      (
      ! peek_css< sequence < zero_plus < exactly <'-' > >, identifier > >(position)
      && lex_css< sequence< negate< digit >, exactly<'-'> > >()
      )

    ) {

      bool right_ws = peek < css_comments >() != NULL;
      operators.push_back({ lexed.to_string() == "+" ? Sass_OP::ADD : Sass_OP::SUB, left_ws, right_ws });
      operands.push_back(parse_operators());
      left_ws = peek < css_comments >() != NULL;
    }

    if (operands.size() == 0) return lhs;
    Expression_Obj ex = fold_operands(lhs, operands, operators);
    state.offset = pstate - state + pstate.offset;
    ex->pstate(state);
    return ex;
  }

  // parse addition and subtraction operations
  Expression_Obj Parser::parse_operators()
  {
    NESTING_GUARD(nestings);
    advanceToNextToken();
    ParserState state(pstate);
    Expression_Obj factor = parse_factor();
    // if it's a singleton, return it (don't wrap it)
    std::vector<Expression_Obj> operands; // factors
    std::vector<Operand> operators; // ops
    // lex operations to apply to lhs
    const char* left_ws = peek < css_comments >();
    while (lex_css< class_char< static_ops > >()) {
      const char* right_ws = peek < css_comments >();
      switch(*lexed.begin) {
        case '*': operators.push_back({ Sass_OP::MUL, left_ws != 0, right_ws != 0 }); break;
        case '/': operators.push_back({ Sass_OP::DIV, left_ws != 0, right_ws != 0 }); break;
        case '%': operators.push_back({ Sass_OP::MOD, left_ws != 0, right_ws != 0 }); break;
        default: throw std::runtime_error("unknown static op parsed");
      }
      operands.push_back(parse_factor());
      left_ws = peek < css_comments >();
    }
    // operands and operators to binary expression
    Expression_Obj ex = fold_operands(factor, operands, operators);
    state.offset = pstate - state + pstate.offset;
    ex->pstate(state);
    return ex;
  }
  // EO parse_operators


  // called from parse_operators
  // called from parse_value_schema
  Expression_Obj Parser::parse_factor()
  {
    NESTING_GUARD(nestings);
    lex < css_comments >(false);
    if (lex_css< exactly<'('> >()) {
      // parse_map may return a list
      Expression_Obj value = parse_map();
      // lex the expected closing parenthesis
      if (!lex_css< exactly<')'> >()) error("unclosed parenthesis");
      // expression can be evaluated
      return value;
    }
    else if (lex_css< exactly<'['> >()) {
      // explicit bracketed
      Expression_Obj value = parse_bracket_list();
      // lex the expected closing square bracket
      if (!lex_css< exactly<']'> >()) error("unclosed squared bracket");
      return value;
    }
    // string may be interpolated
    // if (lex< quoted_string >()) {
    //   return &parse_string();
    // }
    else if (peek< ie_property >()) {
      return parse_ie_property();
    }
    else if (peek< ie_keyword_arg >()) {
      return parse_ie_keyword_arg();
    }
    else if (peek< sequence < calc_fn_call, exactly <'('> > >()) {
      return parse_calc_function();
    }
    else if (lex < functional_schema >()) {
      return parse_function_call_schema();
    }
    else if (lex< identifier_schema >()) {
      String_Obj string = parse_identifier_schema();
      if (String_Schema_Ptr schema = Cast<String_Schema>(string)) {
        if (lex < exactly < '(' > >()) {
          schema->append(parse_list());
          lex < exactly < ')' > >();
        }
      }
      return string;
    }
    else if (peek< sequence< uri_prefix, W, real_uri_value > >()) {
      return parse_url_function_string();
    }
    else if (peek< re_functional >()) {
      return parse_function_call();
    }
    else if (lex< exactly<'+'> >()) {
      Unary_Expression_Ptr ex = SASS_MEMORY_NEW(Unary_Expression, pstate, Unary_Expression::PLUS, parse_factor());
      if (ex && ex->operand()) ex->is_delayed(ex->operand()->is_delayed());
      return ex;
    }
    else if (lex< exactly<'-'> >()) {
      Unary_Expression_Ptr ex = SASS_MEMORY_NEW(Unary_Expression, pstate, Unary_Expression::MINUS, parse_factor());
      if (ex && ex->operand()) ex->is_delayed(ex->operand()->is_delayed());
      return ex;
    }
    else if (lex< exactly<'/'> >()) {
      Unary_Expression_Ptr ex = SASS_MEMORY_NEW(Unary_Expression, pstate, Unary_Expression::SLASH, parse_factor());
      if (ex && ex->operand()) ex->is_delayed(ex->operand()->is_delayed());
      return ex;
    }
    else if (lex< sequence< kwd_not > >()) {
      Unary_Expression_Ptr ex = SASS_MEMORY_NEW(Unary_Expression, pstate, Unary_Expression::NOT, parse_factor());
      if (ex && ex->operand()) ex->is_delayed(ex->operand()->is_delayed());
      return ex;
    }
    // this whole branch is never hit via spec tests
    else if (peek < sequence < one_plus < alternatives < css_whitespace, exactly<'-'>, exactly<'+'> > >, number > >()) {
      if (parse_number_prefix()) return parse_value(); // prefix is positive
      Unary_Expression_Ptr ex = SASS_MEMORY_NEW(Unary_Expression, pstate, Unary_Expression::MINUS, parse_value());
      if (ex->operand()) ex->is_delayed(ex->operand()->is_delayed());
      return ex;
    }
    else {
      return parse_value();
    }
  }

  bool number_has_zero(const std::string& parsed)
  {
    size_t L = parsed.length();
    return !( (L > 0 && parsed.substr(0, 1) == ".") ||
              (L > 1 && parsed.substr(0, 2) == "0.") ||
              (L > 1 && parsed.substr(0, 2) == "-.")  ||
              (L > 2 && parsed.substr(0, 3) == "-0.") );
  }

  Number_Ptr Parser::lexed_number(const ParserState& pstate, const std::string& parsed)
  {
    Number_Ptr nr = SASS_MEMORY_NEW(Number,
                                    pstate,
                                    sass_strtod(parsed.c_str()),
                                    "",
                                    number_has_zero(parsed));
    nr->is_interpolant(false);
    nr->is_delayed(true);
    return nr;
  }

  Number_Ptr Parser::lexed_percentage(const ParserState& pstate, const std::string& parsed)
  {
    Number_Ptr nr = SASS_MEMORY_NEW(Number,
                                    pstate,
                                    sass_strtod(parsed.c_str()),
                                    "%",
                                    true);
    nr->is_interpolant(false);
    nr->is_delayed(true);
    return nr;
  }

  Number_Ptr Parser::lexed_dimension(const ParserState& pstate, const std::string& parsed)
  {
    size_t L = parsed.length();
    size_t num_pos = parsed.find_first_not_of(" \n\r\t");
    if (num_pos == std::string::npos) num_pos = L;
    size_t unit_pos = parsed.find_first_not_of("-+0123456789.", num_pos);
    if (parsed[unit_pos] == 'e' && is_number(parsed[unit_pos+1]) ) {
      unit_pos = parsed.find_first_not_of("-+0123456789.", ++ unit_pos);
    }
    if (unit_pos == std::string::npos) unit_pos = L;
    const std::string& num = parsed.substr(num_pos, unit_pos - num_pos);
    Number_Ptr nr = SASS_MEMORY_NEW(Number,
                                    pstate,
                                    sass_strtod(num.c_str()),
                                    Token(number(parsed.c_str())),
                                    number_has_zero(parsed));
    nr->is_interpolant(false);
    nr->is_delayed(true);
    return nr;
  }

  Value_Ptr Parser::lexed_hex_color(const ParserState& pstate, const std::string& parsed)
  {
    Color_RGBA_Ptr color = NULL;
    if (parsed[0] != '#') {
      return SASS_MEMORY_NEW(String_Quoted, pstate, parsed);
    }
    // chop off the '#'
    std::string hext(parsed.substr(1));
    if (parsed.length() == 4) {
      std::string r(2, parsed[1]);
      std::string g(2, parsed[2]);
      std::string b(2, parsed[3]);
      color = SASS_MEMORY_NEW(Color_RGBA,
                               pstate,
                               static_cast<double>(strtol(r.c_str(), NULL, 16)),
                               static_cast<double>(strtol(g.c_str(), NULL, 16)),
                               static_cast<double>(strtol(b.c_str(), NULL, 16)),
                               1, // alpha channel
                               parsed);
    }
    else if (parsed.length() == 5) {
      std::string r(2, parsed[1]);
      std::string g(2, parsed[2]);
      std::string b(2, parsed[3]);
      std::string a(2, parsed[4]);
      color = SASS_MEMORY_NEW(Color_RGBA,
                               pstate,
                               static_cast<double>(strtol(r.c_str(), NULL, 16)),
                               static_cast<double>(strtol(g.c_str(), NULL, 16)),
                               static_cast<double>(strtol(b.c_str(), NULL, 16)),
                               static_cast<double>(strtol(a.c_str(), NULL, 16)) / 255,
                               parsed);
    }
    else if (parsed.length() == 7) {
      std::string r(parsed.substr(1,2));
      std::string g(parsed.substr(3,2));
      std::string b(parsed.substr(5,2));
      color = SASS_MEMORY_NEW(Color_RGBA,
                               pstate,
                               static_cast<double>(strtol(r.c_str(), NULL, 16)),
                               static_cast<double>(strtol(g.c_str(), NULL, 16)),
                               static_cast<double>(strtol(b.c_str(), NULL, 16)),
                               1, // alpha channel
                               parsed);
    }
    else if (parsed.length() == 9) {
      std::string r(parsed.substr(1,2));
      std::string g(parsed.substr(3,2));
      std::string b(parsed.substr(5,2));
      std::string a(parsed.substr(7,2));
      color = SASS_MEMORY_NEW(Color_RGBA,
                               pstate,
                               static_cast<double>(strtol(r.c_str(), NULL, 16)),
                               static_cast<double>(strtol(g.c_str(), NULL, 16)),
                               static_cast<double>(strtol(b.c_str(), NULL, 16)),
                               static_cast<double>(strtol(a.c_str(), NULL, 16)) / 255,
                               parsed);
    }
    color->is_interpolant(false);
    color->is_delayed(false);
    return color;
  }

  Value_Ptr Parser::color_or_string(const std::string& lexed) const
  {
    if (auto color = name_to_color(lexed)) {
      auto c = SASS_MEMORY_NEW(Color_RGBA, color);
      c->is_delayed(true);
      c->pstate(pstate);
      c->disp(lexed);
      return c;
    } else {
      return SASS_MEMORY_NEW(String_Constant, pstate, lexed);
    }
  }

  // parse one value for a list
  Expression_Obj Parser::parse_value()
  {
    lex< css_comments >(false);
    if (lex< ampersand >())
    {
      if (match< ampersand >()) {
        warning("In Sass, \"&&\" means two copies of the parent selector. You probably want to use \"and\" instead.", pstate);
      }
      return SASS_MEMORY_NEW(Parent_Reference, pstate); }

    if (lex< kwd_important >())
    { return SASS_MEMORY_NEW(String_Constant, pstate, "!important"); }

    // parse `10%4px` into separated items and not a schema
    if (lex< sequence < percentage, lookahead < number > > >())
    { return lexed_percentage(lexed); }

    if (lex< sequence < number, lookahead< sequence < op, number > > > >())
    { return lexed_number(lexed); }

    // string may be interpolated
    if (lex< sequence < quoted_string, lookahead < exactly <'-'> > > >())
    { return parse_string(); }

    if (const char* stop = peek< value_schema >())
    { return parse_value_schema(stop); }

    // string may be interpolated
    if (lex< quoted_string >())
    { return parse_string(); }

    if (lex< kwd_true >())
    { return SASS_MEMORY_NEW(Boolean, pstate, true); }

    if (lex< kwd_false >())
    { return SASS_MEMORY_NEW(Boolean, pstate, false); }

    if (lex< kwd_null >())
    { return SASS_MEMORY_NEW(Null, pstate); }

    if (lex< identifier >()) {
      return color_or_string(lexed);
    }

    if (lex< percentage >())
    { return lexed_percentage(lexed); }

    // match hex number first because 0x000 looks like a number followed by an identifier
    if (lex< sequence < alternatives< hex, hex0 >, negate < exactly<'-'> > > >())
    { return lexed_hex_color(lexed); }

    if (lex< hexa >())
    { return lexed_hex_color(lexed); }

    if (lex< sequence < exactly <'#'>, identifier > >())
    { return SASS_MEMORY_NEW(String_Quoted, pstate, lexed); }

    // also handle the 10em- foo special case
    // alternatives < exactly < '.' >, .. > -- `1.5em-.75em` is split into a list, not a binary expression
    if (lex< sequence< dimension, optional< sequence< exactly<'-'>, lookahead< alternatives < space > > > > > >())
    { return lexed_dimension(lexed); }

    if (lex< sequence< static_component, one_plus< strict_identifier > > >())
    { return SASS_MEMORY_NEW(String_Constant, pstate, lexed); }

    if (lex< number >())
    { return lexed_number(lexed); }

    if (lex< variable >())
    { return SASS_MEMORY_NEW(Variable, pstate, Util::normalize_underscores(lexed)); }

    css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");

    // unreachable statement
    return {};
  }

  // this parses interpolation inside other strings
  // means the result should later be quoted again
  String_Obj Parser::parse_interpolated_chunk(Token chunk, bool constant, bool css)
  {
    const char* i = chunk.begin;
    // see if there any interpolants
    const char* p = constant ? find_first_in_interval< exactly<hash_lbrace> >(i, chunk.end) :
                    find_first_in_interval< exactly<hash_lbrace>, block_comment >(i, chunk.end);

    if (!p) {
      String_Quoted_Ptr str_quoted = SASS_MEMORY_NEW(String_Quoted, pstate, std::string(i, chunk.end), 0, false, false, true, css);
      if (!constant && str_quoted->quote_mark()) str_quoted->quote_mark('*');
      return str_quoted;
    }

    String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate, 0, css);
    schema->is_interpolant(true);
    while (i < chunk.end) {
      p = constant ? find_first_in_interval< exactly<hash_lbrace> >(i, chunk.end) :
          find_first_in_interval< exactly<hash_lbrace>, block_comment >(i, chunk.end);
      if (p) {
        if (i < p) {
          // accumulate the preceding segment if it's nonempty
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string(i, p), css));
        }
        // we need to skip anything inside strings
        // create a new target in parser/prelexer
        if (peek < sequence < optional_spaces, exactly<rbrace> > >(p+2)) { position = p+2;
          css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
        }
        const char* j = skip_over_scopes< exactly<hash_lbrace>, exactly<rbrace> >(p + 2, chunk.end); // find the closing brace
        if (j) { --j;
          // parse the interpolant and accumulate it
          Expression_Obj interp_node = Parser::from_token(Token(p+2, j), ctx, traces, pstate, source).parse_list();
          interp_node->is_interpolant(true);
          schema->append(interp_node);
          i = j;
        }
        else {
          // throw an error if the interpolant is unterminated
          error("unterminated interpolant inside string constant " + chunk.to_string());
        }
      }
      else { // no interpolants left; add the last segment if nonempty
        // check if we need quotes here (was not sure after merge)
        if (i < chunk.end) schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string(i, chunk.end), css));
        break;
      }
      ++ i;
    }

    return schema.detach();
  }

  String_Schema_Obj Parser::parse_css_variable_value(bool top_level)
  {
    String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate);
    String_Schema_Obj tok;
    if (!(tok = parse_css_variable_value_token(top_level))) {
      return {};
    }

    schema->concat(tok);
    while ((tok = parse_css_variable_value_token(top_level))) {
      schema->concat(tok);
    }

    return schema.detach();
  }

  String_Schema_Obj Parser::parse_css_variable_value_token(bool top_level)
  {
    String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate);
    if (
      (top_level && lex< css_variable_top_level_value >(false)) ||
      (!top_level && lex< css_variable_value >(false))
    ) {
      Token str(lexed);
      schema->append(SASS_MEMORY_NEW(String_Constant, pstate, str));
    }
    else if (Expression_Obj tok = lex_interpolation()) {
      if (String_Schema_Ptr s = Cast<String_Schema>(tok)) {
        schema->concat(s);
      } else {
        schema->append(tok);
      }
    }
    else if (lex< quoted_string >()) {
      Expression_Obj tok = parse_string();
      if (String_Schema_Ptr s = Cast<String_Schema>(tok)) {
        schema->concat(s);
      } else {
        schema->append(tok);
      }
    }
    else {
      if (peek< alternatives< exactly<'('>, exactly<'['>, exactly<'{'> > >()) {
        if (lex< exactly<'('> >()) {
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string("(")));
          if (String_Schema_Obj tok = parse_css_variable_value(false)) schema->concat(tok);
          if (!lex< exactly<')'> >()) css_error("Invalid CSS", " after ", ": expected \")\", was ");
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string(")")));
        }
        else if (lex< exactly<'['> >()) {
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string("[")));
          if (String_Schema_Obj tok = parse_css_variable_value(false)) schema->concat(tok);
          if (!lex< exactly<']'> >()) css_error("Invalid CSS", " after ", ": expected \"]\", was ");
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string("]")));
        }
        else if (lex< exactly<'{'> >()) {
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string("{")));
          if (String_Schema_Obj tok = parse_css_variable_value(false)) schema->concat(tok);
          if (!lex< exactly<'}'> >()) css_error("Invalid CSS", " after ", ": expected \"}\", was ");
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string("}")));
        }
      }
    }

    return schema->length() > 0 ? schema.detach() : NULL;
  }

  Value_Obj Parser::parse_static_value()
  {
    lex< static_value >();
    Token str(lexed);
    // static values always have trailing white-
    // space and end delimiter (\s*[;]$) included
    --pstate.offset.column;
    --after_token.column;
    --str.end;
    --position;

    return color_or_string(str.time_wspace());;
  }

  String_Obj Parser::parse_string()
  {
    return parse_interpolated_chunk(Token(lexed));
  }

  String_Obj Parser::parse_ie_property()
  {
    lex< ie_property >();
    Token str(lexed);
    const char* i = str.begin;
    // see if there any interpolants
    const char* p = find_first_in_interval< exactly<hash_lbrace>, block_comment >(str.begin, str.end);
    if (!p) {
      return SASS_MEMORY_NEW(String_Quoted, pstate, std::string(str.begin, str.end));
    }

    String_Schema_Ptr schema = SASS_MEMORY_NEW(String_Schema, pstate);
    while (i < str.end) {
      p = find_first_in_interval< exactly<hash_lbrace>, block_comment >(i, str.end);
      if (p) {
        if (i < p) {
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string(i, p))); // accumulate the preceding segment if it's nonempty
        }
        if (peek < sequence < optional_spaces, exactly<rbrace> > >(p+2)) { position = p+2;
          css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
        }
        const char* j = skip_over_scopes< exactly<hash_lbrace>, exactly<rbrace> >(p+2, str.end); // find the closing brace
        if (j) {
          // parse the interpolant and accumulate it
          Expression_Obj interp_node = Parser::from_token(Token(p+2, j), ctx, traces, pstate, source).parse_list();
          interp_node->is_interpolant(true);
          schema->append(interp_node);
          i = j;
        }
        else {
          // throw an error if the interpolant is unterminated
          error("unterminated interpolant inside IE function " + str.to_string());
        }
      }
      else { // no interpolants left; add the last segment if nonempty
        if (i < str.end) {
          schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string(i, str.end)));
        }
        break;
      }
    }
    return schema;
  }

  String_Obj Parser::parse_ie_keyword_arg()
  {
    String_Schema_Obj kwd_arg = SASS_MEMORY_NEW(String_Schema, pstate, 3);
    if (lex< variable >()) {
      kwd_arg->append(SASS_MEMORY_NEW(Variable, pstate, Util::normalize_underscores(lexed)));
    } else {
      lex< alternatives< identifier_schema, identifier > >();
      kwd_arg->append(SASS_MEMORY_NEW(String_Constant, pstate, lexed));
    }
    lex< exactly<'='> >();
    kwd_arg->append(SASS_MEMORY_NEW(String_Constant, pstate, lexed));
    if (peek< variable >()) kwd_arg->append(parse_list());
    else if (lex< number >()) {
      std::string parsed(lexed);
      Util::normalize_decimals(parsed);
      kwd_arg->append(lexed_number(parsed));
    }
    else if (peek < ie_keyword_arg_value >()) { kwd_arg->append(parse_list()); }
    return kwd_arg;
  }

  String_Schema_Obj Parser::parse_value_schema(const char* stop)
  {
    // initialize the string schema object to add tokens
    String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate);

    if (peek<exactly<'}'>>()) {
      css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
    }

    const char* e;
    const char* ee = end;
    end = stop;
    size_t num_items = 0;
    bool need_space = false;
    while (position < stop) {
      // parse space between tokens
      if (lex< spaces >() && num_items) {
        need_space = true;
      }
      if (need_space) {
        need_space = false;
        // schema->append(SASS_MEMORY_NEW(String_Constant, pstate, " "));
      }
      if ((e = peek< re_functional >()) && e < stop) {
        schema->append(parse_function_call());
      }
      // lex an interpolant /#{...}/
      else if (lex< exactly < hash_lbrace > >()) {
        // Try to lex static expression first
        if (peek< exactly< rbrace > >()) {
          css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
        }
        Expression_Obj ex;
        if (lex< re_static_expression >()) {
          ex = SASS_MEMORY_NEW(String_Constant, pstate, lexed);
        } else {
          ex = parse_list(true);
        }
        ex->is_interpolant(true);
        schema->append(ex);
        if (!lex < exactly < rbrace > >()) {
          css_error("Invalid CSS", " after ", ": expected \"}\", was ");
        }
      }
      // lex some string constants or other valid token
      // Note: [-+] chars are left over from i.e. `#{3}+3`
      else if (lex< alternatives < exactly<'%'>, exactly < '-' >, exactly < '+' > > >()) {
        schema->append(SASS_MEMORY_NEW(String_Constant, pstate, lexed));
      }
      // lex a quoted string
      else if (lex< quoted_string >()) {
        // need_space = true;
        // if (schema->length()) schema->append(SASS_MEMORY_NEW(String_Constant, pstate, " "));
        // else need_space = true;
        schema->append(parse_string());
        if ((*position == '"' || *position == '\'') || peek < alternatives < alpha > >()) {
          // need_space = true;
        }
        if (peek < exactly < '-' > >()) break;
      }
      else if (lex< identifier >()) {
        schema->append(SASS_MEMORY_NEW(String_Constant, pstate, lexed));
        if ((*position == '"' || *position == '\'') || peek < alternatives < alpha > >()) {
           // need_space = true;
        }
      }
      // lex (normalized) variable
      else if (lex< variable >()) {
        std::string name(Util::normalize_underscores(lexed));
        schema->append(SASS_MEMORY_NEW(Variable, pstate, name));
      }
      // lex percentage value
      else if (lex< percentage >()) {
        schema->append(lexed_percentage(lexed));
      }
      // lex dimension value
      else if (lex< dimension >()) {
        schema->append(lexed_dimension(lexed));
      }
      // lex number value
      else if (lex< number >()) {
        schema->append(lexed_number(lexed));
      }
      // lex hex color value
      else if (lex< sequence < hex, negate < exactly < '-' > > > >()) {
        schema->append(lexed_hex_color(lexed));
      }
      else if (lex< sequence < exactly <'#'>, identifier > >()) {
        schema->append(SASS_MEMORY_NEW(String_Quoted, pstate, lexed));
      }
      // lex a value in parentheses
      else if (peek< parenthese_scope >()) {
        schema->append(parse_factor());
      }
      else {
        break;
      }
      ++num_items;
    }
    if (position != stop) {
      schema->append(SASS_MEMORY_NEW(String_Constant, pstate, std::string(position, stop)));
      position = stop;
    }
    end = ee;
    return schema;
  }

  // this parses interpolation outside other strings
  // means the result must not be quoted again later
  String_Obj Parser::parse_identifier_schema()
  {
    Token id(lexed);
    const char* i = id.begin;
    // see if there any interpolants
    const char* p = find_first_in_interval< exactly<hash_lbrace>, block_comment >(id.begin, id.end);
    if (!p) {
      return SASS_MEMORY_NEW(String_Constant, pstate, std::string(id.begin, id.end));
    }

    String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate);
    while (i < id.end) {
      p = find_first_in_interval< exactly<hash_lbrace>, block_comment >(i, id.end);
      if (p) {
        if (i < p) {
          // accumulate the preceding segment if it's nonempty
          const char* o = position; position = i;
          schema->append(parse_value_schema(p));
          position = o;
        }
        // we need to skip anything inside strings
        // create a new target in parser/prelexer
        if (peek < sequence < optional_spaces, exactly<rbrace> > >(p+2)) { position = p;
          css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ");
        }
        const char* j = skip_over_scopes< exactly<hash_lbrace>, exactly<rbrace> >(p+2, id.end); // find the closing brace
        if (j) {
          // parse the interpolant and accumulate it
          Expression_Obj interp_node = Parser::from_token(Token(p+2, j), ctx, traces, pstate, source).parse_list(DELAYED);
          interp_node->is_interpolant(true);
          schema->append(interp_node);
          // schema->has_interpolants(true);
          i = j;
        }
        else {
          // throw an error if the interpolant is unterminated
          error("unterminated interpolant inside interpolated identifier " + id.to_string());
        }
      }
      else { // no interpolants left; add the last segment if nonempty
        if (i < end) {
          const char* o = position; position = i;
          schema->append(parse_value_schema(id.end));
          position = o;
        }
        break;
      }
    }
    return schema ? schema.detach() : 0;
  }

  // calc functions should preserve arguments
  Function_Call_Obj Parser::parse_calc_function()
  {
    lex< identifier >();
    std::string name(lexed);
    ParserState call_pos = pstate;
    lex< exactly<'('> >();
    ParserState arg_pos = pstate;
    const char* arg_beg = position;
    parse_list();
    const char* arg_end = position;
    lex< skip_over_scopes <
          exactly < '(' >,
          exactly < ')' >
        > >();

    Argument_Obj arg = SASS_MEMORY_NEW(Argument, arg_pos, parse_interpolated_chunk(Token(arg_beg, arg_end)));
    Arguments_Obj args = SASS_MEMORY_NEW(Arguments, arg_pos);
    args->append(arg);
    return SASS_MEMORY_NEW(Function_Call, call_pos, name, args);
  }

  String_Obj Parser::parse_url_function_string()
  {
    std::string prefix("");
    if (lex< uri_prefix >()) {
      prefix = std::string(lexed);
    }

    lex < optional_spaces >();
    String_Obj url_string = parse_url_function_argument();

    std::string suffix("");
    if (lex< real_uri_suffix >()) {
      suffix = std::string(lexed);
    }

    std::string uri("");
    if (url_string) {
      uri = url_string->to_string({ NESTED, 5 });
    }

    if (String_Schema_Ptr schema = Cast<String_Schema>(url_string)) {
      String_Schema_Obj res = SASS_MEMORY_NEW(String_Schema, pstate);
      res->append(SASS_MEMORY_NEW(String_Constant, pstate, prefix));
      res->append(schema);
      res->append(SASS_MEMORY_NEW(String_Constant, pstate, suffix));
      return res;
    } else {
      std::string res = prefix + uri + suffix;
      return SASS_MEMORY_NEW(String_Constant, pstate, res);
    }
  }

  String_Obj Parser::parse_url_function_argument()
  {
    const char* p = position;

    std::string uri("");
    if (lex< real_uri_value >(false)) {
      uri = lexed.to_string();
    }

    if (peek< exactly< hash_lbrace > >()) {
      const char* pp = position;
      // TODO: error checking for unclosed interpolants
      while (pp && peek< exactly< hash_lbrace > >(pp)) {
        pp = sequence< interpolant, real_uri_value >(pp);
      }
      if (!pp) return {};
      position = pp;
      return parse_interpolated_chunk(Token(p, position));
    }
    else if (uri != "") {
      std::string res = Util::rtrim(uri);
      return SASS_MEMORY_NEW(String_Constant, pstate, res);
    }

    return {};
  }

  Function_Call_Obj Parser::parse_function_call()
  {
    lex< identifier >();
    std::string name(lexed);

    if (Util::normalize_underscores(name) == "content-exists" && stack.back() != Scope::Mixin)
    { error("Cannot call content-exists() except within a mixin."); }

    ParserState call_pos = pstate;
    Arguments_Obj args = parse_arguments();
    return SASS_MEMORY_NEW(Function_Call, call_pos, name, args);
  }

  Function_Call_Obj Parser::parse_function_call_schema()
  {
    String_Obj name = parse_identifier_schema();
    ParserState source_position_of_call = pstate;
    Arguments_Obj args = parse_arguments();

    return SASS_MEMORY_NEW(Function_Call, source_position_of_call, name, args);
  }

  Content_Obj Parser::parse_content_directive()
  {
    ParserState call_pos = pstate;
    Arguments_Obj args = parse_arguments();

    return SASS_MEMORY_NEW(Content, call_pos, args);
  }

  If_Obj Parser::parse_if_directive(bool else_if)
  {
    stack.push_back(Scope::Control);
    ParserState if_source_position = pstate;
    bool root = block_stack.back()->is_root();
    Expression_Obj predicate = parse_list();
    Block_Obj block = parse_block(root);
    Block_Obj alternative;

    // only throw away comment if we parse a case
    // we want all other comments to be parsed
    if (lex_css< elseif_directive >()) {
      alternative = SASS_MEMORY_NEW(Block, pstate);
      alternative->append(parse_if_directive(true));
    }
    else if (lex_css< kwd_else_directive >()) {
      alternative = parse_block(root);
    }
    stack.pop_back();
    return SASS_MEMORY_NEW(If, if_source_position, predicate, block, alternative);
  }

  For_Obj Parser::parse_for_directive()
  {
    stack.push_back(Scope::Control);
    ParserState for_source_position = pstate;
    bool root = block_stack.back()->is_root();
    lex_variable();
    std::string var(Util::normalize_underscores(lexed));
    if (!lex< kwd_from >()) error("expected 'from' keyword in @for directive");
    Expression_Obj lower_bound = parse_expression();
    bool inclusive = false;
    if (lex< kwd_through >()) inclusive = true;
    else if (lex< kwd_to >()) inclusive = false;
    else                  error("expected 'through' or 'to' keyword in @for directive");
    Expression_Obj upper_bound = parse_expression();
    Block_Obj body = parse_block(root);
    stack.pop_back();
    return SASS_MEMORY_NEW(For, for_source_position, var, lower_bound, upper_bound, body, inclusive);
  }

  // helper to parse a var token
  Token Parser::lex_variable()
  {
    // peek for dollar sign first
    if (!peek< exactly <'$'> >()) {
      css_error("Invalid CSS", " after ", ": expected \"$\", was ");
    }
    // we expect a simple identifier as the call name
    if (!lex< sequence < exactly <'$'>, identifier > >()) {
      lex< exactly <'$'> >(); // move pstate and position up
      css_error("Invalid CSS", " after ", ": expected identifier, was ");
    }
    // return object
    return token;
  }
  // helper to parse identifier
  Token Parser::lex_identifier()
  {
    // we expect a simple identifier as the call name
    if (!lex< identifier >()) { // ToDo: pstate wrong?
      css_error("Invalid CSS", " after ", ": expected identifier, was ");
    }
    // return object
    return token;
  }

  Each_Obj Parser::parse_each_directive()
  {
    stack.push_back(Scope::Control);
    ParserState each_source_position = pstate;
    bool root = block_stack.back()->is_root();
    std::vector<std::string> vars;
    lex_variable();
    vars.push_back(Util::normalize_underscores(lexed));
    while (lex< exactly<','> >()) {
      if (!lex< variable >()) error("@each directive requires an iteration variable");
      vars.push_back(Util::normalize_underscores(lexed));
    }
    if (!lex< kwd_in >()) error("expected 'in' keyword in @each directive");
    Expression_Obj list = parse_list();
    Block_Obj body = parse_block(root);
    stack.pop_back();
    return SASS_MEMORY_NEW(Each, each_source_position, vars, list, body);
  }

  // called after parsing `kwd_while_directive`
  While_Obj Parser::parse_while_directive()
  {
    stack.push_back(Scope::Control);
    bool root = block_stack.back()->is_root();
    // create the initial while call object
    While_Obj call = SASS_MEMORY_NEW(While, pstate, {}, {});
    // parse mandatory predicate
    Expression_Obj predicate = parse_list();
    List_Obj l = Cast<List>(predicate);
    if (!predicate || (l && !l->length())) {
      css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was ", false);
    }
    call->predicate(predicate);
    // parse mandatory block
    call->block(parse_block(root));
    // return ast node
    stack.pop_back();
    // return ast node
    return call.detach();
  }

  // EO parse_while_directive
  Media_Block_Obj Parser::parse_media_block()
  {
    stack.push_back(Scope::Media);
    Media_Block_Obj media_block = SASS_MEMORY_NEW(Media_Block, pstate, {}, {});

    media_block->media_queries(parse_media_queries());

    Media_Block_Obj prev_media_block = last_media_block;
    last_media_block = media_block;
    media_block->block(parse_css_block());
    last_media_block = prev_media_block;
    stack.pop_back();
    return media_block.detach();
  }

  List_Obj Parser::parse_media_queries()
  {
    advanceToNextToken();
    List_Obj queries = SASS_MEMORY_NEW(List, pstate, 0, SASS_COMMA);
    if (!peek_css < exactly <'{'> >()) queries->append(parse_media_query());
    while (lex_css < exactly <','> >()) queries->append(parse_media_query());
    queries->update_pstate(pstate);
    return queries.detach();
  }

  // Expression_Ptr Parser::parse_media_query()
  Media_Query_Obj Parser::parse_media_query()
  {
    advanceToNextToken();
    Media_Query_Obj media_query = SASS_MEMORY_NEW(Media_Query, pstate);
    if (lex < kwd_not >()) { media_query->is_negated(true); lex < css_comments >(false); }
    else if (lex < kwd_only >()) { media_query->is_restricted(true); lex < css_comments >(false); }

    if (lex < identifier_schema >())  media_query->media_type(parse_identifier_schema());
    else if (lex < identifier >())    media_query->media_type(parse_interpolated_chunk(lexed));
    else                             media_query->append(parse_media_expression());

    while (lex_css < kwd_and >()) media_query->append(parse_media_expression());
    if (lex < identifier_schema >()) {
      String_Schema_Ptr schema = SASS_MEMORY_NEW(String_Schema, pstate);
      schema->append(media_query->media_type());
      schema->append(SASS_MEMORY_NEW(String_Constant, pstate, " "));
      schema->append(parse_identifier_schema());
      media_query->media_type(schema);
    }
    while (lex_css < kwd_and >()) media_query->append(parse_media_expression());

    media_query->update_pstate(pstate);

    return media_query;
  }

  Media_Query_Expression_Obj Parser::parse_media_expression()
  {
    if (lex < identifier_schema >()) {
      String_Obj ss = parse_identifier_schema();
      return SASS_MEMORY_NEW(Media_Query_Expression, pstate, ss, {}, true);
    }
    if (!lex_css< exactly<'('> >()) {
      error("media query expression must begin with '('");
    }
    Expression_Obj feature;
    if (peek_css< exactly<')'> >()) {
      error("media feature required in media query expression");
    }
    feature = parse_expression();
    Expression_Obj expression;
    if (lex_css< exactly<':'> >()) {
      expression = parse_list(DELAYED);
    }
    if (!lex_css< exactly<')'> >()) {
      error("unclosed parenthesis in media query expression");
    }
    return SASS_MEMORY_NEW(Media_Query_Expression, feature->pstate(), feature, expression);
  }

  // lexed after `kwd_supports_directive`
  // these are very similar to media blocks
  Supports_Block_Obj Parser::parse_supports_directive()
  {
    Supports_Condition_Obj cond = parse_supports_condition();
    if (!cond) {
      css_error("Invalid CSS", " after ", ": expected @supports condition (e.g. (display: flexbox)), was ", false);
    }
    // create the ast node object for the support queries
    Supports_Block_Obj query = SASS_MEMORY_NEW(Supports_Block, pstate, cond);
    // additional block is mandatory
    // parse inner block
    query->block(parse_block());
    // return ast node
    return query;
  }

  // parse one query operation
  // may encounter nested queries
  Supports_Condition_Obj Parser::parse_supports_condition()
  {
    lex < css_whitespace >();
    Supports_Condition_Obj cond;
    if ((cond = parse_supports_negation())) return cond;
    if ((cond = parse_supports_operator())) return cond;
    if ((cond = parse_supports_interpolation())) return cond;
    return cond;
  }

  Supports_Condition_Obj Parser::parse_supports_negation()
  {
    if (!lex < kwd_not >()) return {};
    Supports_Condition_Obj cond = parse_supports_condition_in_parens();
    return SASS_MEMORY_NEW(Supports_Negation, pstate, cond);
  }

  Supports_Condition_Obj Parser::parse_supports_operator()
  {
    Supports_Condition_Obj cond = parse_supports_condition_in_parens();
    if (cond.isNull()) return {};

    while (true) {
      Supports_Operator::Operand op = Supports_Operator::OR;
      if (lex < kwd_and >()) { op = Supports_Operator::AND; }
      else if(!lex < kwd_or >()) { break; }

      lex < css_whitespace >();
      Supports_Condition_Obj right = parse_supports_condition_in_parens();

      // Supports_Condition_Ptr cc = SASS_MEMORY_NEW(Supports_Condition, *static_cast<Supports_Condition_Ptr>(cond));
      cond = SASS_MEMORY_NEW(Supports_Operator, pstate, cond, right, op);
    }
    return cond;
  }

  Supports_Condition_Obj Parser::parse_supports_interpolation()
  {
    if (!lex < interpolant >()) return {};

    String_Obj interp = parse_interpolated_chunk(lexed);
    if (!interp) return {};

    return SASS_MEMORY_NEW(Supports_Interpolation, pstate, interp);
  }

  // TODO: This needs some major work. Although feature conditions
  // look like declarations their semantics differ significantly
  Supports_Condition_Obj Parser::parse_supports_declaration()
  {
    Supports_Condition_Ptr cond;
    // parse something declaration like
    Expression_Obj feature = parse_expression();
    Expression_Obj expression;
    if (lex_css< exactly<':'> >()) {
      expression = parse_list(DELAYED);
    }
    if (!feature || !expression) error("@supports condition expected declaration");
    cond = SASS_MEMORY_NEW(Supports_Declaration,
                     feature->pstate(),
                     feature,
                     expression);
    // ToDo: maybe we need an additional error condition?
    return cond;
  }

  Supports_Condition_Obj Parser::parse_supports_condition_in_parens()
  {
    Supports_Condition_Obj interp = parse_supports_interpolation();
    if (interp != 0) return interp;

    if (!lex < exactly <'('> >()) return {};
    lex < css_whitespace >();

    Supports_Condition_Obj cond = parse_supports_condition();
    if (cond != 0) {
      if (!lex < exactly <')'> >()) error("unclosed parenthesis in @supports declaration");
    } else {
      cond = parse_supports_declaration();
      if (!lex < exactly <')'> >()) error("unclosed parenthesis in @supports declaration");
    }
    lex < css_whitespace >();
    return cond;
  }

  At_Root_Block_Obj Parser::parse_at_root_block()
  {
    stack.push_back(Scope::AtRoot);
    ParserState at_source_position = pstate;
    Block_Obj body;
    At_Root_Query_Obj expr;
    Lookahead lookahead_result;
    if (lex_css< exactly<'('> >()) {
      expr = parse_at_root_query();
    }
    if (peek_css < exactly<'{'> >()) {
      lex <optional_spaces>();
      body = parse_block(true);
    }
    else if ((lookahead_result = lookahead_for_selector(position)).found) {
      Ruleset_Obj r = parse_ruleset(lookahead_result);
      body = SASS_MEMORY_NEW(Block, r->pstate(), 1, true);
      body->append(r);
    }
    At_Root_Block_Obj at_root = SASS_MEMORY_NEW(At_Root_Block, at_source_position, body);
    if (!expr.isNull()) at_root->expression(expr);
    stack.pop_back();
    return at_root;
  }

  At_Root_Query_Obj Parser::parse_at_root_query()
  {
    if (peek< exactly<')'> >()) error("at-root feature required in at-root expression");

    if (!peek< alternatives< kwd_with_directive, kwd_without_directive > >()) {
      css_error("Invalid CSS", " after ", ": expected \"with\" or \"without\", was ");
    }

    Expression_Obj feature = parse_list();
    if (!lex_css< exactly<':'> >()) error("style declaration must contain a value");
    Expression_Obj expression = parse_list();
    List_Obj value = SASS_MEMORY_NEW(List, feature->pstate(), 1);

    if (expression->concrete_type() == Expression::LIST) {
        value = Cast<List>(expression);
    }
    else value->append(expression);

    At_Root_Query_Obj cond = SASS_MEMORY_NEW(At_Root_Query,
                                          value->pstate(),
                                          feature,
                                          value);
    if (!lex_css< exactly<')'> >()) error("unclosed parenthesis in @at-root expression");
    return cond;
  }

  Directive_Obj Parser::parse_special_directive()
  {
    std::string kwd(lexed);

    if (lexed == "@else") error("Invalid CSS: @else must come after @if");

    // this whole branch is never hit via spec tests

    Directive_Ptr at_rule = SASS_MEMORY_NEW(Directive, pstate, kwd);
    Lookahead lookahead = lookahead_for_include(position);
    if (lookahead.found && !lookahead.has_interpolants) {
      at_rule->selector(parse_selector_list(false));
    }

    lex < css_comments >(false);

    if (lex < static_property >()) {
      at_rule->value(parse_interpolated_chunk(Token(lexed)));
    } else if (!(peek < alternatives < exactly<'{'>, exactly<'}'>, exactly<';'> > >())) {
      at_rule->value(parse_list());
    }

    lex < css_comments >(false);

    if (peek< exactly<'{'> >()) {
      at_rule->block(parse_block());
    }

    return at_rule;
  }

  // this whole branch is never hit via spec tests
  Directive_Obj Parser::parse_prefixed_directive()
  {
    std::string kwd(lexed);

    if (lexed == "@else") error("Invalid CSS: @else must come after @if");

    Directive_Obj at_rule = SASS_MEMORY_NEW(Directive, pstate, kwd);
    Lookahead lookahead = lookahead_for_include(position);
    if (lookahead.found && !lookahead.has_interpolants) {
      at_rule->selector(parse_selector_list(false));
    }

    lex < css_comments >(false);

    if (lex < static_property >()) {
      at_rule->value(parse_interpolated_chunk(Token(lexed)));
    } else if (!(peek < alternatives < exactly<'{'>, exactly<'}'>, exactly<';'> > >())) {
      at_rule->value(parse_list());
    }

    lex < css_comments >(false);

    if (peek< exactly<'{'> >()) {
      at_rule->block(parse_block());
    }

    return at_rule;
  }


  Directive_Obj Parser::parse_directive()
  {
    Directive_Obj directive = SASS_MEMORY_NEW(Directive, pstate, lexed);
    String_Schema_Obj val = parse_almost_any_value();
    // strip left and right if they are of type string
    directive->value(val);
    if (peek< exactly<'{'> >()) {
      directive->block(parse_block());
    }
    return directive;
  }

  Expression_Obj Parser::lex_interpolation()
  {
    if (lex < interpolant >(true) != NULL) {
      return parse_interpolated_chunk(lexed, true);
    }
    return {};
  }

  Expression_Obj Parser::lex_interp_uri()
  {
    // create a string schema by lexing optional interpolations
    return lex_interp< re_string_uri_open, re_string_uri_close >();
  }

  Expression_Obj Parser::lex_interp_string()
  {
    Expression_Obj rv;
    if ((rv = lex_interp< re_string_double_open, re_string_double_close >())) return rv;
    if ((rv = lex_interp< re_string_single_open, re_string_single_close >())) return rv;
    return rv;
  }

  Expression_Obj Parser::lex_almost_any_value_chars()
  {
    const char* match =
    lex <
      one_plus <
        alternatives <
          sequence <
            exactly <'\\'>,
            any_char
          >,
          sequence <
            negate <
              sequence <
                exactly < url_kwd >,
                exactly <'('>
              >
            >,
            neg_class_char <
              almost_any_value_class
            >
          >,
          sequence <
            exactly <'/'>,
            negate <
              alternatives <
                exactly <'/'>,
                exactly <'*'>
              >
            >
          >,
          sequence <
            exactly <'\\'>,
            exactly <'#'>,
            negate <
              exactly <'{'>
            >
          >,
          sequence <
            exactly <'!'>,
            negate <
              alpha
            >
          >
        >
      >
    >(false);
    if (match) {
      return SASS_MEMORY_NEW(String_Constant, pstate, lexed);
    }
    return {};
  }

  Expression_Obj Parser::lex_almost_any_value_token()
  {
    Expression_Obj rv;
    if (*position == 0) return {};
    if ((rv = lex_almost_any_value_chars())) return rv;
    // if ((rv = lex_block_comment())) return rv;
    // if ((rv = lex_single_line_comment())) return rv;
    if ((rv = lex_interp_string())) return rv;
    if ((rv = lex_interp_uri())) return rv;
    if ((rv = lex_interpolation())) return rv;
     if (lex< alternatives< hex, hex0 > >())
    { return lexed_hex_color(lexed); }
   return rv;
  }

  String_Schema_Obj Parser::parse_almost_any_value()
  {

    String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate);
    if (*position == 0) return {};
    lex < spaces >(false);
    Expression_Obj token = lex_almost_any_value_token();
    if (!token) return {};
    schema->append(token);
    if (*position == 0) {
      schema->rtrim();
      return schema.detach();
    }

    while ((token = lex_almost_any_value_token())) {
      schema->append(token);
    }

    lex < css_whitespace >();

    schema->rtrim();

    return schema.detach();
  }

  Warning_Obj Parser::parse_warning()
  {
    if (stack.back() != Scope::Root &&
        stack.back() != Scope::Function &&
        stack.back() != Scope::Mixin &&
        stack.back() != Scope::Control &&
        stack.back() != Scope::Rules) {
      error("Illegal nesting: Only properties may be nested beneath properties.");
    }
    return SASS_MEMORY_NEW(Warning, pstate, parse_list(DELAYED));
  }

  Error_Obj Parser::parse_error()
  {
    if (stack.back() != Scope::Root &&
        stack.back() != Scope::Function &&
        stack.back() != Scope::Mixin &&
        stack.back() != Scope::Control &&
        stack.back() != Scope::Rules) {
      error("Illegal nesting: Only properties may be nested beneath properties.");
    }
    return SASS_MEMORY_NEW(Error, pstate, parse_list(DELAYED));
  }

  Debug_Obj Parser::parse_debug()
  {
    if (stack.back() != Scope::Root &&
        stack.back() != Scope::Function &&
        stack.back() != Scope::Mixin &&
        stack.back() != Scope::Control &&
        stack.back() != Scope::Rules) {
      error("Illegal nesting: Only properties may be nested beneath properties.");
    }
    return SASS_MEMORY_NEW(Debug, pstate, parse_list(DELAYED));
  }

  Return_Obj Parser::parse_return_directive()
  {
    // check that we do not have an empty list (ToDo: check if we got all cases)
    if (peek_css < alternatives < exactly < ';' >, exactly < '}' >, end_of_file > >())
    { css_error("Invalid CSS", " after ", ": expected expression (e.g. 1px, bold), was "); }
    return SASS_MEMORY_NEW(Return, pstate, parse_list());
  }

  Lookahead Parser::lookahead_for_selector(const char* start)
  {
    // init result struct
    Lookahead rv = Lookahead();
    // get start position
    const char* p = start ? start : position;
    // match in one big "regex"
    rv.error = p;
    if (const char* q =
      peek <
        re_selector_list
      >(p)
    ) {
      bool could_be_property = peek< sequence< exactly<'-'>, exactly<'-'> > >(p) != 0;
      bool could_be_escaped = false;
      while (p < q) {
        // did we have interpolations?
        if (*p == '#' && *(p+1) == '{') {
          rv.has_interpolants = true;
          p = q; break;
        }
        // A property that's ambiguous with a nested selector is interpreted as a
        // custom property.
        if (*p == ':' && !could_be_escaped) {
          rv.is_custom_property = could_be_property || p+1 == q || peek< space >(p+1);
        }
        could_be_escaped = *p == '\\';
        ++ p;
      }
      // store anyway  }


      // ToDo: remove
      rv.error = q;
      rv.position = q;
      // check expected opening bracket
      // only after successfull matching
      if (peek < exactly<'{'> >(q)) rv.found = q;
      // else if (peek < end_of_file >(q)) rv.found = q;
      else if (peek < exactly<'('> >(q)) rv.found = q;
      // else if (peek < exactly<';'> >(q)) rv.found = q;
      // else if (peek < exactly<'}'> >(q)) rv.found = q;
      if (rv.found || *p == 0) rv.error = 0;
    }

    rv.parsable = ! rv.has_interpolants;

    // return result
    return rv;

  }
  // EO lookahead_for_selector

  // used in parse_block_nodes and parse_special_directive
  // ToDo: actual usage is still not really clear to me?
  Lookahead Parser::lookahead_for_include(const char* start)
  {
    // we actually just lookahead for a selector
    Lookahead rv = lookahead_for_selector(start);
    // but the "found" rules are different
    if (const char* p = rv.position) {
      // check for additional abort condition
      if (peek < exactly<';'> >(p)) rv.found = p;
      else if (peek < exactly<'}'> >(p)) rv.found = p;
    }
    // return result
    return rv;
  }
  // EO lookahead_for_include

  // look ahead for a token with interpolation in it
  // we mostly use the result if there is an interpolation
  // everything that passes here gets parsed as one schema
  // meaning it will not be parsed as a space separated list
  Lookahead Parser::lookahead_for_value(const char* start)
  {
    // init result struct
    Lookahead rv = Lookahead();
    // get start position
    const char* p = start ? start : position;
    // match in one big "regex"
    if (const char* q =
      peek <
        non_greedy <
          alternatives <
            // consume whitespace
            block_comment, // spaces,
            // main tokens
            sequence <
              interpolant,
              optional <
                quoted_string
              >
            >,
            identifier,
            variable,
            // issue #442
            sequence <
              parenthese_scope,
              interpolant,
              optional <
                quoted_string
              >
            >
          >,
          sequence <
            // optional_spaces,
            alternatives <
              // end_of_file,
              exactly<'{'>,
              exactly<'}'>,
              exactly<';'>
            >
          >
        >
      >(p)
    ) {
      if (p == q) return rv;
      while (p < q) {
        // did we have interpolations?
        if (*p == '#' && *(p+1) == '{') {
          rv.has_interpolants = true;
          p = q; break;
        }
        ++ p;
      }
      // store anyway
      // ToDo: remove
      rv.position = q;
      // check expected opening bracket
      // only after successful matching
      if (peek < exactly<'{'> >(q)) rv.found = q;
      else if (peek < exactly<';'> >(q)) rv.found = q;
      else if (peek < exactly<'}'> >(q)) rv.found = q;
    }

    // return result
    return rv;
  }
  // EO lookahead_for_value

  void Parser::read_bom()
  {
    size_t skip = 0;
    std::string encoding;
    bool utf_8 = false;
    switch ((unsigned char) source[0]) {
    case 0xEF:
      skip = check_bom_chars(source, end, utf_8_bom, 3);
      encoding = "UTF-8";
      utf_8 = true;
      break;
    case 0xFE:
      skip = check_bom_chars(source, end, utf_16_bom_be, 2);
      encoding = "UTF-16 (big endian)";
      break;
    case 0xFF:
      skip = check_bom_chars(source, end, utf_16_bom_le, 2);
      skip += (skip ? check_bom_chars(source, end, utf_32_bom_le, 4) : 0);
      encoding = (skip == 2 ? "UTF-16 (little endian)" : "UTF-32 (little endian)");
      break;
    case 0x00:
      skip = check_bom_chars(source, end, utf_32_bom_be, 4);
      encoding = "UTF-32 (big endian)";
      break;
    case 0x2B:
      skip = check_bom_chars(source, end, utf_7_bom_1, 4)
           | check_bom_chars(source, end, utf_7_bom_2, 4)
           | check_bom_chars(source, end, utf_7_bom_3, 4)
           | check_bom_chars(source, end, utf_7_bom_4, 4)
           | check_bom_chars(source, end, utf_7_bom_5, 5);
      encoding = "UTF-7";
      break;
    case 0xF7:
      skip = check_bom_chars(source, end, utf_1_bom, 3);
      encoding = "UTF-1";
      break;
    case 0xDD:
      skip = check_bom_chars(source, end, utf_ebcdic_bom, 4);
      encoding = "UTF-EBCDIC";
      break;
    case 0x0E:
      skip = check_bom_chars(source, end, scsu_bom, 3);
      encoding = "SCSU";
      break;
    case 0xFB:
      skip = check_bom_chars(source, end, bocu_1_bom, 3);
      encoding = "BOCU-1";
      break;
    case 0x84:
      skip = check_bom_chars(source, end, gb_18030_bom, 4);
      encoding = "GB-18030";
      break;
    default: break;
    }
    if (skip > 0 && !utf_8) error("only UTF-8 documents are currently supported; your document appears to be " + encoding);
    position += skip;
  }

  size_t check_bom_chars(const char* src, const char *end, const unsigned char* bom, size_t len)
  {
    size_t skip = 0;
    if (src + len > end) return 0;
    for (size_t i = 0; i < len; ++i, ++skip) {
      if ((unsigned char) src[i] != bom[i]) return 0;
    }
    return skip;
  }


  Expression_Obj Parser::fold_operands(Expression_Obj base, std::vector<Expression_Obj>& operands, Operand op)
  {
    for (size_t i = 0, S = operands.size(); i < S; ++i) {
      base = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), op, base, operands[i]);
    }
    return base;
  }

  Expression_Obj Parser::fold_operands(Expression_Obj base, std::vector<Expression_Obj>& operands, std::vector<Operand>& ops, size_t i)
  {
    if (String_Schema_Ptr schema = Cast<String_Schema>(base)) {
      // return schema;
      if (schema->has_interpolants()) {
        if (i + 1 < operands.size() && (
             (ops[0].operand == Sass_OP::EQ)
          || (ops[0].operand == Sass_OP::ADD)
          || (ops[0].operand == Sass_OP::DIV)
          || (ops[0].operand == Sass_OP::MUL)
          || (ops[0].operand == Sass_OP::NEQ)
          || (ops[0].operand == Sass_OP::LT)
          || (ops[0].operand == Sass_OP::GT)
          || (ops[0].operand == Sass_OP::LTE)
          || (ops[0].operand == Sass_OP::GTE)
        )) {
          Expression_Obj rhs = fold_operands(operands[i], operands, ops, i + 1);
          rhs = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), ops[0], schema, rhs);
          return rhs;
        }
        // return schema;
      }
    }

    for (size_t S = operands.size(); i < S; ++i) {
      if (String_Schema_Ptr schema = Cast<String_Schema>(operands[i])) {
        if (schema->has_interpolants()) {
          if (i + 1 < S) {
            // this whole branch is never hit via spec tests
            Expression_Obj rhs = fold_operands(operands[i+1], operands, ops, i + 2);
            rhs = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), ops[i], schema, rhs);
            base = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), ops[i], base, rhs);
            return base;
          }
          base = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), ops[i], base, operands[i]);
          return base;
        } else {
          base = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), ops[i], base, operands[i]);
        }
      } else {
        base = SASS_MEMORY_NEW(Binary_Expression, base->pstate(), ops[i], base, operands[i]);
      }
      Binary_Expression_Ptr b = Cast<Binary_Expression>(base.ptr());
      if (b && ops[i].operand == Sass_OP::DIV && b->left()->is_delayed() && b->right()->is_delayed()) {
        base->is_delayed(true);
      }
    }
    // nested binary expression are never to be delayed
    if (Binary_Expression_Ptr b = Cast<Binary_Expression>(base)) {
      if (Cast<Binary_Expression>(b->left())) base->set_delayed(false);
      if (Cast<Binary_Expression>(b->right())) base->set_delayed(false);
    }
    return base;
  }

  void Parser::error(std::string msg, Position pos)
  {
    Position p(pos.line ? pos : before_token);
    ParserState pstate(path, source, p, Offset(0, 0));
    // `pstate.src` may not outlive stack unwind so we must copy it.
    char *src_copy = sass_copy_c_string(pstate.src);
    pstate.src = src_copy;
    traces.push_back(Backtrace(pstate));
    throw Exception::InvalidSass(pstate, traces, msg, src_copy);
  }

  void Parser::error(std::string msg)
  {
    error(msg, pstate);
  }

  // print a css parsing error with actual context information from parsed source
  void Parser::css_error(const std::string& msg, const std::string& prefix, const std::string& middle, const bool trim)
  {
    int max_len = 18;
    const char* end = this->end;
    while (*end != 0) ++ end;
    const char* pos = peek < optional_spaces >();
    if (!pos) pos = position;

    const char* last_pos(pos);
    if (last_pos > source) {
      utf8::prior(last_pos, source);
    }
    // backup position to last significant char
    while (trim && last_pos > source && last_pos < end) {
      if (!Prelexer::is_space(*last_pos)) break;
      utf8::prior(last_pos, source);
    }

    bool ellipsis_left = false;
    const char* pos_left(last_pos);
    const char* end_left(last_pos);

    if (*pos_left) utf8::next(pos_left, end);
    if (*end_left) utf8::next(end_left, end);
    while (pos_left > source) {
      if (utf8::distance(pos_left, end_left) >= max_len) {
        utf8::prior(pos_left, source);
        ellipsis_left = *(pos_left) != '\n' &&
                        *(pos_left) != '\r';
        utf8::next(pos_left, end);
        break;
      }

      const char* prev = pos_left;
      utf8::prior(prev, source);
      if (*prev == '\r') break;
      if (*prev == '\n') break;
      pos_left = prev;
    }
    if (pos_left < source) {
      pos_left = source;
    }

    bool ellipsis_right = false;
    const char* end_right(pos);
    const char* pos_right(pos);
    while (end_right < end) {
      if (utf8::distance(pos_right, end_right) > max_len) {
        ellipsis_left = *(pos_right) != '\n' &&
                        *(pos_right) != '\r';
        break;
      }
      if (*end_right == '\r') break;
      if (*end_right == '\n') break;
      utf8::next(end_right, end);
    }
    // if (*end_right == 0) end_right ++;

    std::string left(pos_left, end_left);
    std::string right(pos_right, end_right);
    size_t left_subpos = left.size() > 15 ? left.size() - 15 : 0;
    size_t right_subpos = right.size() > 15 ? right.size() - 15 : 0;
    if (left_subpos && ellipsis_left) left = ellipsis + left.substr(left_subpos);
    if (right_subpos && ellipsis_right) right = right.substr(right_subpos) + ellipsis;
    // Hotfix when source is null, probably due to interpolation parsing!?
    if (source == NULL || *source == 0) source = pstate.src;
    // now pass new message to the more generic error function
    error(msg + prefix + quote(left) + middle + quote(right));
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Supports_Block::Supports_Block(ParserState pstate, Supports_Condition_Obj condition, Block_Obj block)
  : Has_Block(pstate, block), condition_(condition)
  { statement_type(SUPPORTS); }
  Supports_Block::Supports_Block(const Supports_Block* ptr)
  : Has_Block(ptr), condition_(ptr->condition_)
  { statement_type(SUPPORTS); }
  bool Supports_Block::bubbles() { return true; }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Supports_Condition::Supports_Condition(ParserState pstate)
  : Expression(pstate)
  { }

  Supports_Condition::Supports_Condition(const Supports_Condition* ptr)
  : Expression(ptr)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Supports_Operator::Supports_Operator(ParserState pstate, Supports_Condition_Obj l, Supports_Condition_Obj r, Operand o)
  : Supports_Condition(pstate), left_(l), right_(r), operand_(o)
  { }
  Supports_Operator::Supports_Operator(const Supports_Operator* ptr)
  : Supports_Condition(ptr),
    left_(ptr->left_),
    right_(ptr->right_),
    operand_(ptr->operand_)
  { }

  bool Supports_Operator::needs_parens(Supports_Condition_Obj cond) const
  {
    if (Supports_Operator_Obj op = Cast<Supports_Operator>(cond)) {
      return op->operand() != operand();
    }
    return Cast<Supports_Negation>(cond) != NULL;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Supports_Negation::Supports_Negation(ParserState pstate, Supports_Condition_Obj c)
  : Supports_Condition(pstate), condition_(c)
  { }
  Supports_Negation::Supports_Negation(const Supports_Negation* ptr)
  : Supports_Condition(ptr), condition_(ptr->condition_)
  { }

  bool Supports_Negation::needs_parens(Supports_Condition_Obj cond) const
  {
    return Cast<Supports_Negation>(cond) ||
           Cast<Supports_Operator>(cond);
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Supports_Declaration::Supports_Declaration(ParserState pstate, Expression_Obj f, Expression_Obj v)
  : Supports_Condition(pstate), feature_(f), value_(v)
  { }
  Supports_Declaration::Supports_Declaration(const Supports_Declaration* ptr)
  : Supports_Condition(ptr),
    feature_(ptr->feature_),
    value_(ptr->value_)
  { }

  bool Supports_Declaration::needs_parens(Supports_Condition_Obj cond) const
  {
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Supports_Interpolation::Supports_Interpolation(ParserState pstate, Expression_Obj v)
  : Supports_Condition(pstate), value_(v)
  { }
  Supports_Interpolation::Supports_Interpolation(const Supports_Interpolation* ptr)
  : Supports_Condition(ptr),
    value_(ptr->value_)
  { }

  bool Supports_Interpolation::needs_parens(Supports_Condition_Obj cond) const
  {
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  IMPLEMENT_AST_OPERATORS(Supports_Block);
  IMPLEMENT_AST_OPERATORS(Supports_Condition);
  IMPLEMENT_AST_OPERATORS(Supports_Operator);
  IMPLEMENT_AST_OPERATORS(Supports_Negation);
  IMPLEMENT_AST_OPERATORS(Supports_Declaration);
  IMPLEMENT_AST_OPERATORS(Supports_Interpolation);

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {
  namespace Constants {

    extern const unsigned long MaxCallStack = 1024;

    // https://github.com/sass/libsass/issues/592
    // https://developer.mozilla.org/en-US/docs/Web/CSS/Specificity
    // https://github.com/sass/sass/issues/1495#issuecomment-61189114
    extern const unsigned long Specificity_Star = 0;
    extern const unsigned long Specificity_Universal = 0;
    extern const unsigned long Specificity_Element = 1;
    extern const unsigned long Specificity_Base = 1000;
    extern const unsigned long Specificity_Class = 1000;
    extern const unsigned long Specificity_Attr = 1000;
    extern const unsigned long Specificity_Pseudo = 1000;
    extern const unsigned long Specificity_ID = 1000000;

    extern const int UnificationOrder_Element = 1;
    extern const int UnificationOrder_Id = 2;
    extern const int UnificationOrder_Class = 2;
    extern const int UnificationOrder_Attribute = 3;
    extern const int UnificationOrder_PseudoClass = 4;
    extern const int UnificationOrder_Wrapped = 5;
    extern const int UnificationOrder_PseudoElement = 6;
    extern const int UnificationOrder_Placeholder = 7;

    // sass keywords
    extern const char at_root_kwd[]       = "@at-root";
    extern const char import_kwd[]        = "@import";
    extern const char mixin_kwd[]         = "@mixin";
    extern const char function_kwd[]      = "@function";
    extern const char return_kwd[]        = "@return";
    extern const char include_kwd[]       = "@include";
    extern const char content_kwd[]       = "@content";
    extern const char extend_kwd[]        = "@extend";
    extern const char if_kwd[]            = "@if";
    extern const char else_kwd[]          = "@else";
    extern const char if_after_else_kwd[] = "if";
    extern const char for_kwd[]           = "@for";
    extern const char from_kwd[]          = "from";
    extern const char to_kwd[]            = "to";
    extern const char through_kwd[]       = "through";
    extern const char each_kwd[]          = "@each";
    extern const char in_kwd[]            = "in";
    extern const char while_kwd[]         = "@while";
    extern const char warn_kwd[]          = "@warn";
    extern const char error_kwd[]         = "@error";
    extern const char debug_kwd[]         = "@debug";
    extern const char default_kwd[]       = "default";
    extern const char global_kwd[]        = "global";
    extern const char null_kwd[]          = "null";
    extern const char optional_kwd[]      = "optional";
    extern const char with_kwd[]          = "with";
    extern const char without_kwd[]       = "without";
    extern const char all_kwd[]           = "all";
    extern const char rule_kwd[]          = "rule";

    // css standard units
    extern const char em_kwd[]   = "em";
    extern const char ex_kwd[]   = "ex";
    extern const char px_kwd[]   = "px";
    extern const char cm_kwd[]   = "cm";
    extern const char mm_kwd[]   = "mm";
    extern const char pt_kwd[]   = "pt";
    extern const char pc_kwd[]   = "pc";
    extern const char deg_kwd[]  = "deg";
    extern const char rad_kwd[]  = "rad";
    extern const char grad_kwd[] = "grad";
    extern const char turn_kwd[] = "turn";
    extern const char ms_kwd[]   = "ms";
    extern const char s_kwd[]    = "s";
    extern const char Hz_kwd[]   = "Hz";
    extern const char kHz_kwd[]  = "kHz";

    // vendor prefixes
    extern const char vendor_opera_kwd[]    = "-o-";
    extern const char vendor_webkit_kwd[]   = "-webkit-";
    extern const char vendor_mozilla_kwd[]  = "-moz-";
    extern const char vendor_ms_kwd[]       = "-ms-";
    extern const char vendor_khtml_kwd[]    = "-khtml-";

    // css functions and keywords
    extern const char charset_kwd[]      = "@charset";
    extern const char media_kwd[]        = "@media";
    extern const char supports_kwd[]     = "@supports";
    extern const char keyframes_kwd[]    = "keyframes";
    extern const char only_kwd[]         = "only";
    extern const char rgb_fn_kwd[]       = "rgb(";
    extern const char url_fn_kwd[]       = "url(";
    extern const char url_kwd[]          = "url";
    // extern const char url_prefix_fn_kwd[] = "url-prefix(";
    extern const char important_kwd[]    = "important";
    extern const char pseudo_not_fn_kwd[] = ":not(";
    extern const char even_kwd[]         = "even";
    extern const char odd_kwd[]          = "odd";
    extern const char progid_kwd[]       = "progid";
    extern const char expression_kwd[]   = "expression";
    extern const char calc_fn_kwd[]      = "calc";

    extern const char almost_any_value_class[] = "\"'#!;{}";

    // css selector keywords
    extern const char sel_deep_kwd[] = "/deep/";

    // css attribute-matching operators
    extern const char tilde_equal[]  = "~=";
    extern const char pipe_equal[]   = "|=";
    extern const char caret_equal[]  = "^=";
    extern const char dollar_equal[] = "$=";
    extern const char star_equal[]   = "*=";

    // relational & logical operators and constants
    extern const char and_kwd[]   = "and";
    extern const char or_kwd[]    = "or";
    extern const char not_kwd[]   = "not";
    extern const char gt[]        = ">";
    extern const char gte[]       = ">=";
    extern const char lt[]        = "<";
    extern const char lte[]       = "<=";
    extern const char eq[]        = "==";
    extern const char neq[]       = "!=";
    extern const char true_kwd[]  = "true";
    extern const char false_kwd[] = "false";

    // definition keywords
    extern const char using_kwd[]   = "using";

    // miscellaneous punctuation and delimiters
    extern const char percent_str[]     = "%";
    extern const char empty_str[]       = "";
    extern const char slash_slash[]     = "//";
    extern const char slash_star[]      = "/*";
    extern const char star_slash[]      = "*/";
    extern const char hash_lbrace[]     = "#{";
    extern const char rbrace[]          = "}";
    extern const char rparen[]          = ")";
    extern const char sign_chars[]      = "-+";
    extern const char op_chars[]        = "-+";
    extern const char hyphen[]          = "-";
    extern const char ellipsis[]        = "...";
    // extern const char url_space_chars[] = " \t\r\n\f";
    // type names
    extern const char numeric_name[]    = "numeric value";
    extern const char number_name[]     = "number";
    extern const char percentage_name[] = "percentage";
    extern const char dimension_name[]  = "numeric dimension";
    extern const char string_name[]     = "string";
    extern const char bool_name[]       = "bool";
    extern const char color_name[]      = "color";
    extern const char list_name[]       = "list";
    extern const char map_name[]        = "map";
    extern const char arglist_name[]    = "arglist";

    // constants for uri parsing (RFC 3986 Appendix A.)
    extern const char uri_chars[]  = ":;/?!%&#@|[]{}'`^\"*+-.,_=~";
    extern const char real_uri_chars[]  = "#%&";

    // some specific constant character classes
    // they must be static to be useable by lexer
    extern const char static_ops[]      = "*/%";
    // some character classes for the parser
    extern const char selector_list_delims[] = "){};!";
    extern const char complex_selector_delims[] = ",){};!";
    extern const char selector_combinator_ops[] = "+~>";
    // optional modifiers for alternative compare context
    extern const char attribute_compare_modifiers[] = "~|^$*";
    extern const char selector_lookahead_ops[] = "*&%,()[]";

    // byte order marks
    // (taken from http://en.wikipedia.org/wiki/Byte_order_mark)
    extern const unsigned char utf_8_bom[]      = { 0xEF, 0xBB, 0xBF };
    extern const unsigned char utf_16_bom_be[]  = { 0xFE, 0xFF };
    extern const unsigned char utf_16_bom_le[]  = { 0xFF, 0xFE };
    extern const unsigned char utf_32_bom_be[]  = { 0x00, 0x00, 0xFE, 0xFF };
    extern const unsigned char utf_32_bom_le[]  = { 0xFF, 0xFE, 0x00, 0x00 };
    extern const unsigned char utf_7_bom_1[]    = { 0x2B, 0x2F, 0x76, 0x38 };
    extern const unsigned char utf_7_bom_2[]    = { 0x2B, 0x2F, 0x76, 0x39 };
    extern const unsigned char utf_7_bom_3[]    = { 0x2B, 0x2F, 0x76, 0x2B };
    extern const unsigned char utf_7_bom_4[]    = { 0x2B, 0x2F, 0x76, 0x2F };
    extern const unsigned char utf_7_bom_5[]    = { 0x2B, 0x2F, 0x76, 0x38, 0x2D };
    extern const unsigned char utf_1_bom[]      = { 0xF7, 0x64, 0x4C };
    extern const unsigned char utf_ebcdic_bom[] = { 0xDD, 0x73, 0x66, 0x73 };
    extern const unsigned char scsu_bom[]       = { 0x0E, 0xFE, 0xFF };
    extern const unsigned char bocu_1_bom[]     = { 0xFB, 0xEE, 0x28 };
    extern const unsigned char gb_18030_bom[]   = { 0x84, 0x31, 0x95, 0x33 };

  }
}

namespace Sass {

  Value_Ptr c2ast(union Sass_Value* v, Backtraces traces, ParserState pstate)
  {
    using std::strlen;
    using std::strcpy;
    Value_Ptr e = NULL;
    switch (sass_value_get_tag(v)) {
      case SASS_BOOLEAN: {
        e = SASS_MEMORY_NEW(Boolean, pstate, !!sass_boolean_get_value(v));
      } break;
      case SASS_NUMBER: {
        e = SASS_MEMORY_NEW(Number, pstate, sass_number_get_value(v), sass_number_get_unit(v));
      } break;
      case SASS_COLOR: {
        e = SASS_MEMORY_NEW(Color_RGBA, pstate, sass_color_get_r(v), sass_color_get_g(v), sass_color_get_b(v), sass_color_get_a(v));
      } break;
      case SASS_STRING: {
        if (sass_string_is_quoted(v))
          e = SASS_MEMORY_NEW(String_Quoted, pstate, sass_string_get_value(v));
        else {
          e = SASS_MEMORY_NEW(String_Constant, pstate, sass_string_get_value(v));
        }
      } break;
      case SASS_LIST: {
        List_Ptr l = SASS_MEMORY_NEW(List, pstate, sass_list_get_length(v), sass_list_get_separator(v));
        for (size_t i = 0, L = sass_list_get_length(v); i < L; ++i) {
          l->append(c2ast(sass_list_get_value(v, i), traces, pstate));
        }
        l->is_bracketed(sass_list_get_is_bracketed(v));
        e = l;
      } break;
      case SASS_MAP: {
        Map_Ptr m = SASS_MEMORY_NEW(Map, pstate);
        for (size_t i = 0, L = sass_map_get_length(v); i < L; ++i) {
          *m << std::make_pair(
            c2ast(sass_map_get_key(v, i), traces, pstate),
            c2ast(sass_map_get_value(v, i), traces, pstate));
        }
        e = m;
      } break;
      case SASS_NULL: {
        e = SASS_MEMORY_NEW(Null, pstate);
      } break;
      case SASS_ERROR: {
        error("Error in C function: " + std::string(sass_error_get_message(v)), pstate, traces);
      } break;
      case SASS_WARNING: {
        error("Warning in C function: " + std::string(sass_warning_get_message(v)), pstate, traces);
      } break;
      default: break;
    }
    return e;
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {
  using namespace Constants;
  using namespace File;
  using namespace Sass;

  inline bool sort_importers (const Sass_Importer_Entry& i, const Sass_Importer_Entry& j)
  { return sass_importer_get_priority(i) > sass_importer_get_priority(j); }

  static std::string safe_input(const char* in_path)
  {
    // enforce some safe defaults
    // used to create relative file links
    std::string safe_path(in_path ? in_path : "");
    return safe_path == "" ? "stdin" : safe_path;
  }

  static std::string safe_output(const char* out_path, const std::string& input_path = "")
  {
    std::string safe_path(out_path ? out_path : "");
    // maybe we can extract an output path from input path
    if (safe_path == "" && input_path != "") {
      int lastindex = static_cast<int>(input_path.find_last_of("."));
      return (lastindex > -1 ? input_path.substr(0, lastindex) : input_path) + ".css";
    }
    // enforce some safe defaults
    // used to create relative file links
    return safe_path == "" ? "stdout" : safe_path;
  }

  Context::Context(struct Sass_Context& c_ctx)
  : CWD(File::get_cwd()),
    c_options(c_ctx),
    entry_path(""),
    head_imports(0),
    plugins(),
    emitter(c_options),

    ast_gc(),
    strings(),
    resources(),
    sheets(),
    subset_map(),
    import_stack(),
    callee_stack(),
    traces(),
    c_compiler(NULL),

    c_headers               (std::vector<Sass_Importer_Entry>()),
    c_importers             (std::vector<Sass_Importer_Entry>()),
    c_functions             (std::vector<Sass_Function_Entry>()),

    indent                  (safe_str(c_options.indent, "  ")),
    linefeed                (safe_str(c_options.linefeed, "\n")),

    input_path              (make_canonical_path(safe_input(c_options.input_path))),
    output_path             (make_canonical_path(safe_output(c_options.output_path, input_path))),
    source_map_file         (make_canonical_path(safe_str(c_options.source_map_file, ""))),
    source_map_root         (make_canonical_path(safe_str(c_options.source_map_root, "")))

  {

    // Sass 3.4: The current working directory will no longer be placed onto the Sass load path by default.
    // If you need the current working directory to be available, set SASS_PATH=. in your shell's environment.
    // include_paths.push_back(CWD);

    // collect more paths from different options
    collect_include_paths(c_options.include_path);
    collect_include_paths(c_options.include_paths);
    collect_plugin_paths(c_options.plugin_path);
    collect_plugin_paths(c_options.plugin_paths);

    // load plugins and register custom behaviors
    for(auto plug : plugin_paths) plugins.load_plugins(plug);
    for(auto fn : plugins.get_headers()) c_headers.push_back(fn);
    for(auto fn : plugins.get_importers()) c_importers.push_back(fn);
    for(auto fn : plugins.get_functions()) c_functions.push_back(fn);

    // sort the items by priority (lowest first)
    sort (c_headers.begin(), c_headers.end(), sort_importers);
    sort (c_importers.begin(), c_importers.end(), sort_importers);

    emitter.set_filename(abs2rel(output_path, source_map_file, CWD));

  }

  void Context::add_c_function(Sass_Function_Entry function)
  {
    c_functions.push_back(function);
  }
  void Context::add_c_header(Sass_Importer_Entry header)
  {
    c_headers.push_back(header);
    // need to sort the array afterwards (no big deal)
    sort (c_headers.begin(), c_headers.end(), sort_importers);
  }
  void Context::add_c_importer(Sass_Importer_Entry importer)
  {
    c_importers.push_back(importer);
    // need to sort the array afterwards (no big deal)
    sort (c_importers.begin(), c_importers.end(), sort_importers);
  }

  Context::~Context()
  {
    // resources were allocated by malloc
    for (size_t i = 0; i < resources.size(); ++i) {
      free(resources[i].contents);
      free(resources[i].srcmap);
    }
    // free all strings we kept alive during compiler execution
    for (size_t n = 0; n < strings.size(); ++n) free(strings[n]);
    // everything that gets put into sources will be freed by us
    // this shouldn't have anything in it anyway!?
    for (size_t m = 0; m < import_stack.size(); ++m) {
      sass_import_take_source(import_stack[m]);
      sass_import_take_srcmap(import_stack[m]);
      sass_delete_import(import_stack[m]);
    }
    // clear inner structures (vectors) and input source
    resources.clear(); import_stack.clear();
    subset_map.clear(), sheets.clear();
  }

  Data_Context::~Data_Context()
  {
    // --> this will be freed by resources
    // make sure we free the source even if not processed!
    // if (resources.size() == 0 && source_c_str) free(source_c_str);
    // if (resources.size() == 0 && srcmap_c_str) free(srcmap_c_str);
    // source_c_str = 0; srcmap_c_str = 0;
  }

  File_Context::~File_Context()
  {
  }

  void Context::collect_include_paths(const char* paths_str)
  {
    if (paths_str) {
      const char* beg = paths_str;
      const char* end = Prelexer::find_first<PATH_SEP>(beg);

      while (end) {
        std::string path(beg, end - beg);
        if (!path.empty()) {
          if (*path.rbegin() != '/') path += '/';
          include_paths.push_back(path);
        }
        beg = end + 1;
        end = Prelexer::find_first<PATH_SEP>(beg);
      }

      std::string path(beg);
      if (!path.empty()) {
        if (*path.rbegin() != '/') path += '/';
        include_paths.push_back(path);
      }
    }
  }

  void Context::collect_include_paths(string_list* paths_array)
  {
    while (paths_array)
    {
      collect_include_paths(paths_array->string);
      paths_array = paths_array->next;
    }
  }

  void Context::collect_plugin_paths(const char* paths_str)
  {
    if (paths_str) {
      const char* beg = paths_str;
      const char* end = Prelexer::find_first<PATH_SEP>(beg);

      while (end) {
        std::string path(beg, end - beg);
        if (!path.empty()) {
          if (*path.rbegin() != '/') path += '/';
          plugin_paths.push_back(path);
        }
        beg = end + 1;
        end = Prelexer::find_first<PATH_SEP>(beg);
      }

      std::string path(beg);
      if (!path.empty()) {
        if (*path.rbegin() != '/') path += '/';
        plugin_paths.push_back(path);
      }
    }
  }

  void Context::collect_plugin_paths(string_list* paths_array)
  {
    while (paths_array)
    {
      collect_plugin_paths(paths_array->string);
      paths_array = paths_array->next;
    }
  }

  // resolve the imp_path in base_path or include_paths
  // looks for alternatives and returns a list from one directory
  std::vector<Include> Context::find_includes(const Importer& import)
  {
    // make sure we resolve against an absolute path
    std::string base_path(rel2abs(import.base_path));
    // first try to resolve the load path relative to the base path
    std::vector<Include> vec(resolve_includes(base_path, import.imp_path));
    // then search in every include path (but only if nothing found yet)
    for (size_t i = 0, S = include_paths.size(); vec.size() == 0 && i < S; ++i)
    {
      // call resolve_includes and individual base path and append all results
      std::vector<Include> resolved(resolve_includes(include_paths[i], import.imp_path));
      if (resolved.size()) vec.insert(vec.end(), resolved.begin(), resolved.end());
    }
    // return vector
    return vec;
  }

  // register include with resolved path and its content
  // memory of the resources will be freed by us on exit
  void Context::register_resource(const Include& inc, const Resource& res)
  {

    // do not parse same resource twice
    // maybe raise an error in this case
    // if (sheets.count(inc.abs_path)) {
    //   free(res.contents); free(res.srcmap);
    //   throw std::runtime_error("duplicate resource registered");
    //   return;
    // }

    // get index for this resource
    size_t idx = resources.size();

    // tell emitter about new resource
    emitter.add_source_index(idx);

    // put resources under our control
    // the memory will be freed later
    resources.push_back(res);

    // add a relative link to the working directory
    included_files.push_back(inc.abs_path);
    // add a relative link  to the source map output file
    srcmap_links.push_back(abs2rel(inc.abs_path, source_map_file, CWD));

    // get pointer to the loaded content
    Sass_Import_Entry import = sass_make_import(
      inc.imp_path.c_str(),
      inc.abs_path.c_str(),
      res.contents,
      res.srcmap
    );
    // add the entry to the stack
    import_stack.push_back(import);

    // get pointer to the loaded content
    const char* contents = resources[idx].contents;
    // keep a copy of the path around (for parserstates)
    // ToDo: we clean it, but still not very elegant!?
    strings.push_back(sass_copy_c_string(inc.abs_path.c_str()));
    // create the initial parser state from resource
    ParserState pstate(strings.back(), contents, idx);

    // check existing import stack for possible recursion
    for (size_t i = 0; i < import_stack.size() - 2; ++i) {
      auto parent = import_stack[i];
      if (std::strcmp(parent->abs_path, import->abs_path) == 0) {
        std::string cwd(File::get_cwd());
        // make path relative to the current directory
        std::string stack("An @import loop has been found:");
        for (size_t n = 1; n < i + 2; ++n) {
          stack += "\n    " + std::string(File::abs2rel(import_stack[n]->abs_path, cwd, cwd)) +
            " imports " + std::string(File::abs2rel(import_stack[n+1]->abs_path, cwd, cwd));
        }
        // implement error throw directly until we
        // decided how to handle full stack traces
        throw Exception::InvalidSyntax(pstate, traces, stack);
        // error(stack, prstate ? *prstate : pstate, import_stack);
      }
    }

    // create a parser instance from the given c_str buffer
    Parser p(Parser::from_c_str(contents, *this, traces, pstate));
    // do not yet dispose these buffers
    sass_import_take_source(import);
    sass_import_take_srcmap(import);
    // then parse the root block
    Block_Obj root = p.parse();
    // delete memory of current stack frame
    sass_delete_import(import_stack.back());
    // remove current stack frame
    import_stack.pop_back();
    // create key/value pair for ast node
    std::pair<const std::string, StyleSheet>
      ast_pair(inc.abs_path, { res, root });
    // register resulting resource
    sheets.insert(ast_pair);
  }

  // register include with resolved path and its content
  // memory of the resources will be freed by us on exit
  void Context::register_resource(const Include& inc, const Resource& res, ParserState& prstate)
  {
    traces.push_back(Backtrace(prstate));
    register_resource(inc, res);
    traces.pop_back();
  }

  // Add a new import to the context (called from `import_url`)
  Include Context::load_import(const Importer& imp, ParserState pstate)
  {

    // search for valid imports (ie. partials) on the filesystem
    // this may return more than one valid result (ambiguous imp_path)
    const std::vector<Include> resolved(find_includes(imp));

    // error nicely on ambiguous imp_path
    if (resolved.size() > 1) {
      std::stringstream msg_stream;
      msg_stream << "It's not clear which file to import for ";
      msg_stream << "'@import \"" << imp.imp_path << "\"'." << "\n";
      msg_stream << "Candidates:" << "\n";
      for (size_t i = 0, L = resolved.size(); i < L; ++i)
      { msg_stream << "  " << resolved[i].imp_path << "\n"; }
      msg_stream << "Please delete or rename all but one of these files." << "\n";
      error(msg_stream.str(), pstate, traces);
    }

    // process the resolved entry
    else if (resolved.size() == 1) {
      bool use_cache = c_importers.size() == 0;
      // use cache for the resource loading
      if (use_cache && sheets.count(resolved[0].abs_path)) return resolved[0];
      // try to read the content of the resolved file entry
      // the memory buffer returned must be freed by us!
      if (char* contents = read_file(resolved[0].abs_path)) {
        // register the newly resolved file resource
        register_resource(resolved[0], { contents, 0 }, pstate);
        // return resolved entry
        return resolved[0];
      }
    }

    // nothing found
    return { imp, "" };

  }

  void Context::import_url (Import_Ptr imp, std::string load_path, const std::string& ctx_path) {

    ParserState pstate(imp->pstate());
    std::string imp_path(unquote(load_path));
    std::string protocol("file");

    using namespace Prelexer;
    if (const char* proto = sequence< identifier, exactly<':'>, exactly<'/'>, exactly<'/'> >(imp_path.c_str())) {

      protocol = std::string(imp_path.c_str(), proto - 3);
      // if (protocol.compare("file") && true) { }
    }

    // add urls (protocol other than file) and urls without procotol to `urls` member
    // ToDo: if ctx_path is already a file resource, we should not add it here?
    if (imp->import_queries() || protocol != "file" || imp_path.substr(0, 2) == "//") {
      imp->urls().push_back(SASS_MEMORY_NEW(String_Quoted, imp->pstate(), load_path));
    }
    else if (imp_path.length() > 4 && imp_path.substr(imp_path.length() - 4, 4) == ".css") {
      String_Constant_Ptr loc = SASS_MEMORY_NEW(String_Constant, pstate, unquote(load_path));
      Argument_Obj loc_arg = SASS_MEMORY_NEW(Argument, pstate, loc);
      Arguments_Obj loc_args = SASS_MEMORY_NEW(Arguments, pstate);
      loc_args->append(loc_arg);
      Function_Call_Ptr new_url = SASS_MEMORY_NEW(Function_Call, pstate, std::string("url"), loc_args);
      imp->urls().push_back(new_url);
    }
    else {
      const Importer importer(imp_path, ctx_path);
      Include include(load_import(importer, pstate));
      if (include.abs_path.empty()) {
        error("File to import not found or unreadable: " + imp_path + ".", pstate, traces);
      }
      imp->incs().push_back(include);
    }

  }


  // call custom importers on the given (unquoted) load_path and eventually parse the resulting style_sheet
  bool Context::call_loader(const std::string& load_path, const char* ctx_path, ParserState& pstate, Import_Ptr imp, std::vector<Sass_Importer_Entry> importers, bool only_one)
  {
    // unique counter
    size_t count = 0;
    // need one correct import
    bool has_import = false;
    // process all custom importers (or custom headers)
    for (Sass_Importer_Entry& importer_ent : importers) {
      // int priority = sass_importer_get_priority(importer);
      Sass_Importer_Fn fn = sass_importer_get_function(importer_ent);
      // skip importer if it returns NULL
      if (Sass_Import_List includes =
          fn(load_path.c_str(), importer_ent, c_compiler)
      ) {
        // get c pointer copy to iterate over
        Sass_Import_List it_includes = includes;
        while (*it_includes) { ++count;
          // create unique path to use as key
          std::string uniq_path = load_path;
          if (!only_one && count) {
            std::stringstream path_strm;
            path_strm << uniq_path << ":" << count;
            uniq_path = path_strm.str();
          }
          // create the importer struct
          Importer importer(uniq_path, ctx_path);
          // query data from the current include
          Sass_Import_Entry include_ent = *it_includes;
          char* source = sass_import_take_source(include_ent);
          char* srcmap = sass_import_take_srcmap(include_ent);
          size_t line = sass_import_get_error_line(include_ent);
          size_t column = sass_import_get_error_column(include_ent);
          const char *abs_path = sass_import_get_abs_path(include_ent);
          // handle error message passed back from custom importer
          // it may (or may not) override the line and column info
          if (const char* err_message = sass_import_get_error_message(include_ent)) {
            if (source || srcmap) register_resource({ importer, uniq_path }, { source, srcmap }, pstate);
            if (line == std::string::npos && column == std::string::npos) error(err_message, pstate, traces);
            else error(err_message, ParserState(ctx_path, source, Position(line, column)), traces);
          }
          // content for import was set
          else if (source) {
            // resolved abs_path should be set by custom importer
            // use the created uniq_path as fallback (maybe enforce)
            std::string path_key(abs_path ? abs_path : uniq_path);
            // create the importer struct
            Include include(importer, path_key);
            // attach information to AST node
            imp->incs().push_back(include);
            // register the resource buffers
            register_resource(include, { source, srcmap }, pstate);
          }
          // only a path was retuned
          // try to load it like normal
          else if(abs_path) {
            // checks some urls to preserve
            // `http://`, `https://` and `//`
            // or dispatchs to `import_file`
            // which will check for a `.css` extension
            // or resolves the file on the filesystem
            // added and resolved via `add_file`
            // finally stores everything on `imp`
            import_url(imp, abs_path, ctx_path);
          }
          // move to next
          ++it_includes;
        }
        // deallocate the returned memory
        sass_delete_import_list(includes);
        // set success flag
        has_import = true;
        // break out of loop
        if (only_one) break;
      }
    }
    // return result
    return has_import;
  }

  void register_function(Context&, Signature sig, Native_Function f, Env* env);
  void register_function(Context&, Signature sig, Native_Function f, size_t arity, Env* env);
  void register_overload_stub(Context&, std::string name, Env* env);
  void register_built_in_functions(Context&, Env* env);
  void register_c_functions(Context&, Env* env, Sass_Function_List);
  void register_c_function(Context&, Env* env, Sass_Function_Entry);

  char* Context::render(Block_Obj root)
  {
    // check for valid block
    if (!root) return 0;
    // start the render process
    root->perform(&emitter);
    // finish emitter stream
    emitter.finalize();
    // get the resulting buffer from stream
    OutputBuffer emitted = emitter.get_buffer();
    // should we append a source map url?
    if (!c_options.omit_source_map_url) {
      // generate an embeded source map
      if (c_options.source_map_embed) {
        emitted.buffer += linefeed;
        emitted.buffer += format_embedded_source_map();
      }
      // or just link the generated one
      else if (source_map_file != "") {
        emitted.buffer += linefeed;
        emitted.buffer += format_source_mapping_url(source_map_file);
      }
    }
    // create a copy of the resulting buffer string
    // this must be freed or taken over by implementor
    return sass_copy_c_string(emitted.buffer.c_str());
  }

  void Context::apply_custom_headers(Block_Obj root, const char* ctx_path, ParserState pstate)
  {
    // create a custom import to resolve headers
    Import_Obj imp = SASS_MEMORY_NEW(Import, pstate);
    // dispatch headers which will add custom functions
    // custom headers are added to the import instance
    call_headers(entry_path, ctx_path, pstate, imp);
    // increase head count to skip later
    head_imports += resources.size() - 1;
    // add the statement if we have urls
    if (!imp->urls().empty()) root->append(imp);
    // process all other resources (add Import_Stub nodes)
    for (size_t i = 0, S = imp->incs().size(); i < S; ++i) {
      root->append(SASS_MEMORY_NEW(Import_Stub, pstate, imp->incs()[i]));
    }
  }

  Block_Obj File_Context::parse()
  {

    // check if entry file is given
    if (input_path.empty()) return {};

    // create absolute path from input filename
    // ToDo: this should be resolved via custom importers
    std::string abs_path(rel2abs(input_path, CWD));

    // try to load the entry file
    char* contents = read_file(abs_path);

    // alternatively also look inside each include path folder
    // I think this differs from ruby sass (IMO too late to remove)
    for (size_t i = 0, S = include_paths.size(); contents == 0 && i < S; ++i) {
      // build absolute path for this include path entry
      abs_path = rel2abs(input_path, include_paths[i]);
      // try to load the resulting path
      contents = read_file(abs_path);
    }

    // abort early if no content could be loaded (various reasons)
    if (!contents) throw std::runtime_error("File to read not found or unreadable: " + input_path);

    // store entry path
    entry_path = abs_path;

    // create entry only for import stack
    Sass_Import_Entry import = sass_make_import(
      input_path.c_str(),
      entry_path.c_str(),
      contents,
      0
    );
    // add the entry to the stack
    import_stack.push_back(import);

    // create the source entry for file entry
    register_resource({{ input_path, "." }, abs_path }, { contents, 0 });

    // create root ast tree node
    return compile();

  }

  Block_Obj Data_Context::parse()
  {

    // check if source string is given
    if (!source_c_str) return {};

    // convert indented sass syntax
    if(c_options.is_indented_syntax_src) {
      // call sass2scss to convert the string
      char * converted = sass2scss(source_c_str,
        // preserve the structure as much as possible
        SASS2SCSS_PRETTIFY_1 | SASS2SCSS_KEEP_COMMENT);
      // replace old source_c_str with converted
      free(source_c_str); source_c_str = converted;
    }

    // remember entry path (defaults to stdin for string)
    entry_path = input_path.empty() ? "stdin" : input_path;

    // ToDo: this may be resolved via custom importers
    std::string abs_path(rel2abs(entry_path));
    char* abs_path_c_str = sass_copy_c_string(abs_path.c_str());
    strings.push_back(abs_path_c_str);

    // create entry only for the import stack
    Sass_Import_Entry import = sass_make_import(
      entry_path.c_str(),
      abs_path_c_str,
      source_c_str,
      srcmap_c_str
    );
    // add the entry to the stack
    import_stack.push_back(import);

    // register a synthetic resource (path does not really exist, skip in includes)
    register_resource({{ input_path, "." }, input_path }, { source_c_str, srcmap_c_str });

    // create root ast tree node
    return compile();
  }



  // parse root block from includes
  Block_Obj Context::compile()
  {
    // abort if there is no data
    if (resources.size() == 0) return {};
    // get root block from the first style sheet
    Block_Obj root = sheets.at(entry_path).root;
    // abort on invalid root
    if (root.isNull()) return {};
    Env global; // create root environment
    // register built-in functions on env
    register_built_in_functions(*this, &global);
    // register custom functions (defined via C-API)
    for (size_t i = 0, S = c_functions.size(); i < S; ++i)
    { register_c_function(*this, &global, c_functions[i]); }
    // create initial backtrace entry
    // create crtp visitor objects
    Expand expand(*this, &global);
    Cssize cssize(*this);
    CheckNesting check_nesting;
    // check nesting in all files
    for (auto sheet : sheets) {
      auto styles = sheet.second;
      check_nesting(styles.root);
    }
    // expand and eval the tree
    root = expand(root);
    // check nesting
    check_nesting(root);
    // merge and bubble certain rules
    root = cssize(root);
    // should we extend something?
    if (!subset_map.empty()) {
      // create crtp visitor object
      Extend extend(subset_map);
      extend.setEval(expand.eval);
      // extend tree nodes
      extend(root);
    }

    // clean up by removing empty placeholders
    // ToDo: maybe we can do this somewhere else?
    Remove_Placeholders remove_placeholders;
    root->perform(&remove_placeholders);
    // return processed tree
    return root;
  }
  // EO compile

  std::string Context::format_embedded_source_map()
  {
    std::string map = emitter.render_srcmap(*this);
    std::istringstream is( map );
    std::ostringstream buffer;
    base64::encoder E;
    E.encode(is, buffer);
    std::string url = "data:application/json;base64," + buffer.str();
    url.erase(url.size() - 1);
    return "/*# sourceMappingURL=" + url + " */";
  }

  std::string Context::format_source_mapping_url(const std::string& file)
  {
    std::string url = abs2rel(file, output_path, CWD);
    return "/*# sourceMappingURL=" + url + " */";
  }

  char* Context::render_srcmap()
  {
    if (source_map_file == "") return 0;
    std::string map = emitter.render_srcmap(*this);
    return sass_copy_c_string(map.c_str());
  }


  // for data context we want to start after "stdin"
  // we probably always want to skip the header includes?
  std::vector<std::string> Context::get_included_files(bool skip, size_t headers)
  {
      // create a copy of the vector for manipulations
      std::vector<std::string> includes = included_files;
      if (includes.size() == 0) return includes;
      if (skip) { includes.erase( includes.begin(), includes.begin() + 1 + headers); }
      else { includes.erase( includes.begin() + 1, includes.begin() + 1 + headers); }
      includes.erase( std::unique( includes.begin(), includes.end() ), includes.end() );
      std::sort( includes.begin() + (skip ? 0 : 1), includes.end() );
      return includes;
  }

  void register_function(Context& ctx, Signature sig, Native_Function f, Env* env)
  {
    Definition_Ptr def = make_native_function(sig, f, ctx);
    def->environment(env);
    (*env)[def->name() + "[f]"] = def;
  }

  void register_function(Context& ctx, Signature sig, Native_Function f, size_t arity, Env* env)
  {
    Definition_Ptr def = make_native_function(sig, f, ctx);
    std::stringstream ss;
    ss << def->name() << "[f]" << arity;
    def->environment(env);
    (*env)[ss.str()] = def;
  }

  void register_overload_stub(Context& ctx, std::string name, Env* env)
  {
    Definition_Ptr stub = SASS_MEMORY_NEW(Definition,
                                       ParserState("[built-in function]"),
                                       0,
                                       name,
                                       {},
                                       0,
                                       true);
    (*env)[name + "[f]"] = stub;
  }


  void register_built_in_functions(Context& ctx, Env* env)
  {
    using namespace Functions;
    // RGB Functions
    register_function(ctx, rgb_sig, rgb, env);
    register_overload_stub(ctx, "rgba", env);
    register_function(ctx, rgba_4_sig, rgba_4, 4, env);
    register_function(ctx, rgba_2_sig, rgba_2, 2, env);
    register_function(ctx, red_sig, red, env);
    register_function(ctx, green_sig, green, env);
    register_function(ctx, blue_sig, blue, env);
    register_function(ctx, mix_sig, mix, env);
    // HSL Functions
    register_function(ctx, hsl_sig, hsl, env);
    register_function(ctx, hsla_sig, hsla, env);
    register_function(ctx, hue_sig, hue, env);
    register_function(ctx, saturation_sig, saturation, env);
    register_function(ctx, lightness_sig, lightness, env);
    register_function(ctx, adjust_hue_sig, adjust_hue, env);
    register_function(ctx, lighten_sig, lighten, env);
    register_function(ctx, darken_sig, darken, env);
    register_function(ctx, saturate_sig, saturate, env);
    register_function(ctx, desaturate_sig, desaturate, env);
    register_function(ctx, grayscale_sig, grayscale, env);
    register_function(ctx, complement_sig, complement, env);
    register_function(ctx, invert_sig, invert, env);
    // Opacity Functions
    register_function(ctx, alpha_sig, alpha, env);
    register_function(ctx, opacity_sig, alpha, env);
    register_function(ctx, opacify_sig, opacify, env);
    register_function(ctx, fade_in_sig, opacify, env);
    register_function(ctx, transparentize_sig, transparentize, env);
    register_function(ctx, fade_out_sig, transparentize, env);
    // Other Color Functions
    register_function(ctx, adjust_color_sig, adjust_color, env);
    register_function(ctx, scale_color_sig, scale_color, env);
    register_function(ctx, change_color_sig, change_color, env);
    register_function(ctx, ie_hex_str_sig, ie_hex_str, env);
    // String Functions
    register_function(ctx, unquote_sig, sass_unquote, env);
    register_function(ctx, quote_sig, sass_quote, env);
    register_function(ctx, str_length_sig, str_length, env);
    register_function(ctx, str_insert_sig, str_insert, env);
    register_function(ctx, str_index_sig, str_index, env);
    register_function(ctx, str_slice_sig, str_slice, env);
    register_function(ctx, to_upper_case_sig, to_upper_case, env);
    register_function(ctx, to_lower_case_sig, to_lower_case, env);
    // Number Functions
    register_function(ctx, percentage_sig, percentage, env);
    register_function(ctx, round_sig, round, env);
    register_function(ctx, ceil_sig, ceil, env);
    register_function(ctx, floor_sig, floor, env);
    register_function(ctx, abs_sig, abs, env);
    register_function(ctx, min_sig, min, env);
    register_function(ctx, max_sig, max, env);
    register_function(ctx, random_sig, random, env);
    // List Functions
    register_function(ctx, length_sig, length, env);
    register_function(ctx, nth_sig, nth, env);
    register_function(ctx, set_nth_sig, set_nth, env);
    register_function(ctx, index_sig, index, env);
    register_function(ctx, join_sig, join, env);
    register_function(ctx, append_sig, append, env);
    register_function(ctx, zip_sig, zip, env);
    register_function(ctx, list_separator_sig, list_separator, env);
    register_function(ctx, is_bracketed_sig, is_bracketed, env);
    // Map Functions
    register_function(ctx, map_get_sig, map_get, env);
    register_function(ctx, map_merge_sig, map_merge, env);
    register_function(ctx, map_remove_sig, map_remove, env);
    register_function(ctx, map_keys_sig, map_keys, env);
    register_function(ctx, map_values_sig, map_values, env);
    register_function(ctx, map_has_key_sig, map_has_key, env);
    register_function(ctx, keywords_sig, keywords, env);
    // Introspection Functions
    register_function(ctx, type_of_sig, type_of, env);
    register_function(ctx, unit_sig, unit, env);
    register_function(ctx, unitless_sig, unitless, env);
    register_function(ctx, comparable_sig, comparable, env);
    register_function(ctx, variable_exists_sig, variable_exists, env);
    register_function(ctx, global_variable_exists_sig, global_variable_exists, env);
    register_function(ctx, function_exists_sig, function_exists, env);
    register_function(ctx, mixin_exists_sig, mixin_exists, env);
    register_function(ctx, feature_exists_sig, feature_exists, env);
    register_function(ctx, call_sig, call, env);
    register_function(ctx, content_exists_sig, content_exists, env);
    register_function(ctx, get_function_sig, get_function, env);
    // Boolean Functions
    register_function(ctx, not_sig, sass_not, env);
    register_function(ctx, if_sig, sass_if, env);
    // Misc Functions
    register_function(ctx, inspect_sig, inspect, env);
    register_function(ctx, unique_id_sig, unique_id, env);
    // Selector functions
    register_function(ctx, selector_nest_sig, selector_nest, env);
    register_function(ctx, selector_append_sig, selector_append, env);
    register_function(ctx, selector_extend_sig, selector_extend, env);
    register_function(ctx, selector_replace_sig, selector_replace, env);
    register_function(ctx, selector_unify_sig, selector_unify, env);
    register_function(ctx, is_superselector_sig, is_superselector, env);
    register_function(ctx, simple_selectors_sig, simple_selectors, env);
    register_function(ctx, selector_parse_sig, selector_parse, env);
  }

  void register_c_functions(Context& ctx, Env* env, Sass_Function_List descrs)
  {
    while (descrs && *descrs) {
      register_c_function(ctx, env, *descrs);
      ++descrs;
    }
  }
  void register_c_function(Context& ctx, Env* env, Sass_Function_Entry descr)
  {
    Definition_Ptr def = make_c_function(descr, ctx);
    def->environment(env);
    (*env)[def->name() + "[f]"] = def;
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  std::string Base64VLQ::encode(const int number) const
  {
    std::string encoded = "";

    int vlq = to_vlq_signed(number);

    do {
      int digit = vlq & VLQ_BASE_MASK;
      vlq >>= VLQ_BASE_SHIFT;
      if (vlq > 0) {
        digit |= VLQ_CONTINUATION_BIT;
      }
      encoded += base64_encode(digit);
    } while (vlq > 0);

    return encoded;
  }

  char Base64VLQ::base64_encode(const int number) const
  {
    int index = number;
    if (index < 0) index = 0;
    if (index > 63) index = 63;
    return CHARACTERS[index];
  }

  int Base64VLQ::to_vlq_signed(const int number) const
  {
    return (number < 0) ? ((-number) << 1) + 1 : (number << 1) + 0;
  }

  const char* Base64VLQ::CHARACTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

  const int Base64VLQ::VLQ_BASE_SHIFT = 5;
  const int Base64VLQ::VLQ_BASE = 1 << VLQ_BASE_SHIFT;
  const int Base64VLQ::VLQ_BASE_MASK = VLQ_BASE - 1;
  const int Base64VLQ::VLQ_CONTINUATION_BIT = VLQ_BASE;

}


// #define DEBUG_UNIFY

namespace Sass {

  Compound_Selector_Ptr Compound_Selector::unify_with(Compound_Selector_Ptr rhs)
  {
    #ifdef DEBUG_UNIFY
    const std::string debug_call = "unify(Compound[" + this->to_string() + "], Compound[" + rhs->to_string() + "])";
    std::cerr << debug_call << std::endl;
    #endif

    if (empty()) return rhs;
    Compound_Selector_Obj unified = SASS_MEMORY_COPY(rhs);
    for (const Simple_Selector_Obj& sel : elements()) {
      unified = sel->unify_with(unified);
      if (unified.isNull()) break;
    }

    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " = Compound[" << unified->to_string() << "]" << std::endl;
    #endif
    return unified.detach();
  }

  Compound_Selector_Ptr Simple_Selector::unify_with(Compound_Selector_Ptr rhs)
  {
    #ifdef DEBUG_UNIFY
    const std::string debug_call = "unify(Simple[" + this->to_string() + "], Compound[" + rhs->to_string() + "])";
    std::cerr << debug_call << std::endl;
    #endif

    if (rhs->length() == 1) {
      if (rhs->at(0)->is_universal()) {
        Compound_Selector_Ptr this_compound = SASS_MEMORY_NEW(Compound_Selector, pstate(), 1);
        this_compound->append(SASS_MEMORY_COPY(this));
        Compound_Selector_Ptr unified = rhs->at(0)->unify_with(this_compound);
        if (unified == nullptr || unified != this_compound) delete this_compound;

        #ifdef DEBUG_UNIFY
        std::cerr << "> " << debug_call << " = " << "Compound[" << unified->to_string() << "]" << std::endl;
        #endif
        return unified;
      }
    }
    for (const Simple_Selector_Obj& sel : rhs->elements()) {
      if (*this == *sel) {
        #ifdef DEBUG_UNIFY
        std::cerr << "> " << debug_call << " = " << "Compound[" << rhs->to_string() << "]" << std::endl;
        #endif
        return rhs;
      }
    }
    const int lhs_order = this->unification_order();
    size_t i = rhs->length();
    while (i > 0 && lhs_order < rhs->at(i - 1)->unification_order()) --i;
    rhs->insert(rhs->begin() + i, this);
    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " = " << "Compound[" << rhs->to_string() << "]" << std::endl;
    #endif
    return rhs;
  }

  Simple_Selector_Ptr Type_Selector::unify_with(Simple_Selector_Ptr rhs)
  {
    #ifdef DEBUG_UNIFY
    const std::string debug_call = "unify(Type[" + this->to_string() + "], Simple[" + rhs->to_string() + "])";
    std::cerr << debug_call << std::endl;
    #endif

    bool rhs_ns = false;
    if (!(is_ns_eq(*rhs) || rhs->is_universal_ns())) {
      if (!is_universal_ns()) {
        #ifdef DEBUG_UNIFY
        std::cerr << "> " << debug_call << " = nullptr" << std::endl;
        #endif
        return nullptr;
      }
      rhs_ns = true;
    }
    bool rhs_name = false;
    if (!(name_ == rhs->name() || rhs->is_universal())) {
      if (!(is_universal())) {
        #ifdef DEBUG_UNIFY
        std::cerr << "> " << debug_call << " = nullptr" << std::endl;
        #endif
        return nullptr;
      }
      rhs_name = true;
    }
    if (rhs_ns) {
      ns(rhs->ns());
      has_ns(rhs->has_ns());
    }
    if (rhs_name) name(rhs->name());
    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " = Simple[" << this->to_string() << "]" << std::endl;
    #endif
    return this;
  }

  Compound_Selector_Ptr Type_Selector::unify_with(Compound_Selector_Ptr rhs)
  {
    #ifdef DEBUG_UNIFY
    const std::string debug_call = "unify(Type[" + this->to_string() + "], Compound[" + rhs->to_string() + "])";
    std::cerr << debug_call << std::endl;
    #endif

    if (rhs->empty()) {
      rhs->append(this);
      #ifdef DEBUG_UNIFY
      std::cerr << "> " << debug_call << " = Compound[" << rhs->to_string() << "]" << std::endl;
      #endif
      return rhs;
    }
    Type_Selector_Ptr rhs_0 = Cast<Type_Selector>(rhs->at(0));
    if (rhs_0 != nullptr) {
      Simple_Selector_Ptr unified = unify_with(rhs_0);
      if (unified == nullptr) {
        #ifdef DEBUG_UNIFY
        std::cerr << "> " << debug_call << " = nullptr" << std::endl;
        #endif
        return nullptr;
      }
      rhs->elements()[0] = unified;
    } else if (!is_universal() || (has_ns_ && ns_ != "*")) {
      rhs->insert(rhs->begin(), this);
    }
    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " = Compound[" << rhs->to_string() << "]" << std::endl;
    #endif
    return rhs;
  }

  Compound_Selector_Ptr Class_Selector::unify_with(Compound_Selector_Ptr rhs)
  {
    rhs->has_line_break(has_line_break());
    return Simple_Selector::unify_with(rhs);
  }

  Compound_Selector_Ptr Id_Selector::unify_with(Compound_Selector_Ptr rhs)
  {
    for (const Simple_Selector_Obj& sel : rhs->elements()) {
      if (Id_Selector_Ptr id_sel = Cast<Id_Selector>(sel)) {
        if (id_sel->name() != name()) return nullptr;
      }
    }
    rhs->has_line_break(has_line_break());
    return Simple_Selector::unify_with(rhs);
  }

  Compound_Selector_Ptr Pseudo_Selector::unify_with(Compound_Selector_Ptr rhs)
  {
    if (is_pseudo_element()) {
      for (const Simple_Selector_Obj& sel : rhs->elements()) {
        if (Pseudo_Selector_Ptr pseudo_sel = Cast<Pseudo_Selector>(sel)) {
          if (pseudo_sel->is_pseudo_element() && pseudo_sel->name() != name()) return nullptr;
        }
      }
    }
    return Simple_Selector::unify_with(rhs);
  }

  Selector_List_Ptr Complex_Selector::unify_with(Complex_Selector_Ptr rhs)
  {
    #ifdef DEBUG_UNIFY
    const std::string debug_call = "unify(Complex[" + this->to_string() + "], Complex[" + rhs->to_string() + "])";
    std::cerr << debug_call << std::endl;
    #endif

    // get last tails (on the right side)
    Complex_Selector_Ptr l_last = this->mutable_last();
    Complex_Selector_Ptr r_last = rhs->mutable_last();

    // check valid pointers (assertion)
    SASS_ASSERT(l_last, "lhs is null");
    SASS_ASSERT(r_last, "rhs is null");

    // Not sure about this check, but closest way I could check
    // was to see if this is a ruby 'SimpleSequence' equivalent.
    // It seems to do the job correctly as some specs react to this
    if (l_last->combinator() != Combinator::ANCESTOR_OF) return nullptr;
    if (r_last->combinator() != Combinator::ANCESTOR_OF) return nullptr;

    // get the headers for the last tails
    Compound_Selector_Ptr l_last_head = l_last->head();
    Compound_Selector_Ptr r_last_head = r_last->head();

    // check valid head pointers (assertion)
    SASS_ASSERT(l_last_head, "lhs head is null");
    SASS_ASSERT(r_last_head, "rhs head is null");

    // get the unification of the last compound selectors
    Compound_Selector_Obj unified = r_last_head->unify_with(l_last_head);

    // abort if we could not unify heads
    if (unified == nullptr) return nullptr;

    // move the head
    if (l_last_head->is_universal()) l_last->head({});
    r_last->head(unified);

    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " before weave: lhs=" << this->to_string() << " rhs=" << rhs->to_string() << std::endl;
    #endif

    // create nodes from both selectors
    Node lhsNode = complexSelectorToNode(this);
    Node rhsNode = complexSelectorToNode(rhs);

    // Complex_Selector_Obj fake = unified->to_complex();
    // Node unified_node = complexSelectorToNode(fake);
    // // add to permutate the list?
    // rhsNode.plus(unified_node);

    // do some magic we inherit from node and extend
    Node node = subweave(lhsNode, rhsNode);
    Selector_List_Obj result = SASS_MEMORY_NEW(Selector_List, pstate(), node.collection()->size());
    for (auto &item : *node.collection()) {
      result->append(nodeToComplexSelector(Node::naiveTrim(item)));
    }

    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " = " << result->to_string() << std::endl;
    #endif

    // only return if list has some entries
    return result->length() ? result.detach() : nullptr;
  }

  Selector_List_Ptr Selector_List::unify_with(Selector_List_Ptr rhs) {
    #ifdef DEBUG_UNIFY
    const std::string debug_call = "unify(List[" + this->to_string() + "], List[" + rhs->to_string() + "])";
    std::cerr << debug_call << std::endl;
    #endif

    std::vector<Complex_Selector_Obj> result;
    // Unify all of children with RHS's children, storing the results in `unified_complex_selectors`
    for (Complex_Selector_Obj& seq1 : elements()) {
      for (Complex_Selector_Obj& seq2 : rhs->elements()) {
        Complex_Selector_Obj seq1_copy = SASS_MEMORY_CLONE(seq1);
        Complex_Selector_Obj seq2_copy = SASS_MEMORY_CLONE(seq2);
        Selector_List_Obj unified = seq1_copy->unify_with(seq2_copy);
        if (unified) {
          result.reserve(result.size() + unified->length());
          std::copy(unified->begin(), unified->end(), std::back_inserter(result));
        }
      }
    }

    // Creates the final Selector_List by combining all the complex selectors
    Selector_List_Ptr final_result = SASS_MEMORY_NEW(Selector_List, pstate(), result.size());
    for (Complex_Selector_Obj& sel : result) {
      final_result->append(sel);
    }
    #ifdef DEBUG_UNIFY
    std::cerr << "> " << debug_call << " = " << final_result->to_string() << std::endl;
    #endif
    return final_result;
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  namespace Operators {

    inline double add(double x, double y) { return x + y; }
    inline double sub(double x, double y) { return x - y; }
    inline double mul(double x, double y) { return x * y; }
    inline double div(double x, double y) { return x / y; } // x/0 checked by caller

    inline double mod(double x, double y) { // x/0 checked by caller
      if ((x > 0 && y < 0) || (x < 0 && y > 0)) {
        double ret = std::fmod(x, y);
        return ret ? ret + y : ret;
      } else {
        return std::fmod(x, y);
      }
    }

    typedef double (*bop)(double, double);
    bop ops[Sass_OP::NUM_OPS] = {
      0, 0, // and, or
      0, 0, 0, 0, 0, 0, // eq, neq, gt, gte, lt, lte
      add, sub, mul, div, mod
    };

    /* static function, has no pstate or traces */
    bool eq(Expression_Obj lhs, Expression_Obj rhs)
    {
      // operation is undefined if one is not a number
      if (!lhs || !rhs) throw Exception::UndefinedOperation(lhs, rhs, Sass_OP::EQ);
      // use compare operator from ast node
      return *lhs == *rhs;
    }

    /* static function, throws OperationError, has no pstate or traces */
    bool cmp(Expression_Obj lhs, Expression_Obj rhs, const Sass_OP op)
    {
      // can only compare numbers!?
      Number_Obj l = Cast<Number>(lhs);
      Number_Obj r = Cast<Number>(rhs);
      // operation is undefined if one is not a number
      if (!l || !r) throw Exception::UndefinedOperation(lhs, rhs, op);
      // use compare operator from ast node
      return *l < *r;
    }

    /* static functions, throws OperationError, has no pstate or traces */
    bool lt(Expression_Obj lhs, Expression_Obj rhs) { return cmp(lhs, rhs, Sass_OP::LT); }
    bool neq(Expression_Obj lhs, Expression_Obj rhs) { return eq(lhs, rhs) == false; }
    bool gt(Expression_Obj lhs, Expression_Obj rhs) { return !cmp(lhs, rhs, Sass_OP::GT) && neq(lhs, rhs); }
    bool lte(Expression_Obj lhs, Expression_Obj rhs) { return cmp(lhs, rhs, Sass_OP::LTE) || eq(lhs, rhs); }
    bool gte(Expression_Obj lhs, Expression_Obj rhs) { return !cmp(lhs, rhs, Sass_OP::GTE) || eq(lhs, rhs); }

    /* colour math deprecation warning */
    void op_color_deprecation(enum Sass_OP op, std::string lsh, std::string rhs, const ParserState& pstate)
    {
      deprecated(
        "The operation `" + lsh + " " + sass_op_to_name(op) + " " + rhs +
        "` is deprecated and will be an error in future versions.",
        "Consider using Sass's color functions instead.\n"
        "http://sass-lang.com/documentation/Sass/Script/Functions.html#other_color_functions",
        /*with_column=*/false, pstate);
    }

    /* static function, throws OperationError, has no traces but optional pstate for returned value */
    Value_Ptr op_strings(Sass::Operand operand, Value& lhs, Value& rhs, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed)
    {
      enum Sass_OP op = operand.operand;

      String_Quoted_Ptr lqstr = Cast<String_Quoted>(&lhs);
      String_Quoted_Ptr rqstr = Cast<String_Quoted>(&rhs);

      std::string lstr(lqstr ? lqstr->value() : lhs.to_string(opt));
      std::string rstr(rqstr ? rqstr->value() : rhs.to_string(opt));

      if (Cast<Null>(&lhs)) throw Exception::InvalidNullOperation(&lhs, &rhs, op);
      if (Cast<Null>(&rhs)) throw Exception::InvalidNullOperation(&lhs, &rhs, op);

      std::string sep;
      switch (op) {
        case Sass_OP::ADD: sep = "";   break;
        case Sass_OP::SUB: sep = "-";  break;
        case Sass_OP::DIV: sep = "/";  break;
        case Sass_OP::EQ:  sep = "=="; break;
        case Sass_OP::NEQ: sep = "!="; break;
        case Sass_OP::LT:  sep = "<";  break;
        case Sass_OP::GT:  sep = ">";  break;
        case Sass_OP::LTE: sep = "<="; break;
        case Sass_OP::GTE: sep = ">="; break;
        default:
          throw Exception::UndefinedOperation(&lhs, &rhs, op);
        break;
      }

      if (op == Sass_OP::ADD) {
        // create string that might be quoted on output (but do not unquote what we pass)
        return SASS_MEMORY_NEW(String_Quoted, pstate, lstr + rstr, 0, false, true);
      }

      // add whitespace around operator
      // but only if result is not delayed
      if (sep != "" && delayed == false) {
        if (operand.ws_before) sep = " " + sep;
        if (operand.ws_after) sep = sep + " ";
      }

      if (op == Sass_OP::SUB || op == Sass_OP::DIV) {
        if (lqstr && lqstr->quote_mark()) lstr = quote(lstr);
        if (rqstr && rqstr->quote_mark()) rstr = quote(rstr);
      }

      return SASS_MEMORY_NEW(String_Constant, pstate, lstr + sep + rstr);
    }

    /* ToDo: allow to operate also with hsla colors */
    /* static function, throws OperationError, has no traces but optional pstate for returned value */
    Value_Ptr op_colors(enum Sass_OP op, const Color_RGBA& lhs, const Color_RGBA& rhs, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed)
    {

      if (lhs.a() != rhs.a()) {
        throw Exception::AlphaChannelsNotEqual(&lhs, &rhs, op);
      }
      if ((op == Sass_OP::DIV || op == Sass_OP::MOD) && (!rhs.r() || !rhs.g() || !rhs.b())) {
        throw Exception::ZeroDivisionError(lhs, rhs);
      }

      op_color_deprecation(op, lhs.to_string(), rhs.to_string(), pstate);

      return SASS_MEMORY_NEW(Color_RGBA,
                             pstate,
                             ops[op](lhs.r(), rhs.r()),
                             ops[op](lhs.g(), rhs.g()),
                             ops[op](lhs.b(), rhs.b()),
                             lhs.a());
    }

    /* static function, throws OperationError, has no traces but optional pstate for returned value */
    Value_Ptr op_numbers(enum Sass_OP op, const Number& lhs, const Number& rhs, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed)
    {
      double lval = lhs.value();
      double rval = rhs.value();

      if (op == Sass_OP::MOD && rval == 0) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "NaN");
      }

      if (op == Sass_OP::DIV && rval == 0) {
        std::string result(lval ? "Infinity" : "NaN");
        return SASS_MEMORY_NEW(String_Quoted, pstate, result);
      }

      size_t l_n_units = lhs.numerators.size();
      size_t l_d_units = lhs.numerators.size();
      size_t r_n_units = rhs.denominators.size();
      size_t r_d_units = rhs.denominators.size();
      // optimize out the most common and simplest case
      if (l_n_units == r_n_units && l_d_units == r_d_units) {
        if (l_n_units + l_d_units <= 1 && r_n_units + r_d_units <= 1) {
          if (lhs.numerators == rhs.numerators) {
            if (lhs.denominators == rhs.denominators) {
              Number_Ptr v = SASS_MEMORY_COPY(&lhs);
              v->value(ops[op](lval, rval));
              return v;
            }
          }
        }
      }

      Number_Obj v = SASS_MEMORY_COPY(&lhs);

      if (lhs.is_unitless() && (op == Sass_OP::ADD || op == Sass_OP::SUB || op == Sass_OP::MOD)) {
        v->numerators = rhs.numerators;
        v->denominators = rhs.denominators;
      }

      if (op == Sass_OP::MUL) {
        v->value(ops[op](lval, rval));
        v->numerators.insert(v->numerators.end(),
          rhs.numerators.begin(), rhs.numerators.end()
        );
        v->denominators.insert(v->denominators.end(),
          rhs.denominators.begin(), rhs.denominators.end()
        );
        v->reduce();
      }
      else if (op == Sass_OP::DIV) {
        v->value(ops[op](lval, rval));
        v->numerators.insert(v->numerators.end(),
          rhs.denominators.begin(), rhs.denominators.end()
        );
        v->denominators.insert(v->denominators.end(),
          rhs.numerators.begin(), rhs.numerators.end()
        );
        v->reduce();
      }
      else {
        Number ln(lhs), rn(rhs);
        ln.reduce(); rn.reduce();
        double f(rn.convert_factor(ln));
        v->value(ops[op](lval, rn.value() * f));
      }

      v->pstate(pstate);
      return v.detach();
    }

    /* static function, throws OperationError, has no traces but optional pstate for returned value */
    Value_Ptr op_number_color(enum Sass_OP op, const Number& lhs, const Color_RGBA& rhs, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed)
    {
      double lval = lhs.value();

      switch (op) {
        case Sass_OP::ADD:
        case Sass_OP::MUL: {
          op_color_deprecation(op, lhs.to_string(), rhs.to_string(opt), pstate);
          return SASS_MEMORY_NEW(Color_RGBA,
                                pstate,
                                ops[op](lval, rhs.r()),
                                ops[op](lval, rhs.g()),
                                ops[op](lval, rhs.b()),
                                rhs.a());
        }
        case Sass_OP::SUB:
        case Sass_OP::DIV: {
          std::string color(rhs.to_string(opt));
          op_color_deprecation(op, lhs.to_string(), color, pstate);
          return SASS_MEMORY_NEW(String_Quoted,
                                pstate,
                                lhs.to_string(opt)
                                + sass_op_separator(op)
                                + color);
        }
        default: break;
      }
      throw Exception::UndefinedOperation(&lhs, &rhs, op);
    }

    /* static function, throws OperationError, has no traces but optional pstate for returned value */
    Value_Ptr op_color_number(enum Sass_OP op, const Color_RGBA& lhs, const Number& rhs, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed)
    {
      double rval = rhs.value();

      if ((op == Sass_OP::DIV || op == Sass_OP::DIV) && rval == 0) {
        // comparison of Fixnum with Float failed?
        throw Exception::ZeroDivisionError(lhs, rhs);
      }

      op_color_deprecation(op, lhs.to_string(), rhs.to_string(), pstate);

      return SASS_MEMORY_NEW(Color_RGBA,
                            pstate,
                            ops[op](lhs.r(), rval),
                            ops[op](lhs.g(), rval),
                            ops[op](lhs.b(), rval),
                            lhs.a());
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


/*
 NOTES:

 - The print* functions print to cerr. This allows our testing frameworks (like sass-spec) to ignore the output, which
   is very helpful when debugging. The format of the output is mainly to wrap things in square brackets to match what
   ruby already outputs (to make comparisons easier).

 - For the direct porting effort, we're trying to port method-for-method until we get all the tests passing.
   Where applicable, I've tried to include the ruby code above the function for reference until all our tests pass.
   The ruby code isn't always directly portable, so I've tried to include any modified ruby code that was actually
   used for the porting.

 - DO NOT try to optimize yet. We get a tremendous benefit out of comparing the output of each stage of the extend to the ruby
   output at the same stage. This makes it much easier to determine where problems are. Try to keep as close to
   the ruby code as you can until we have all the sass-spec tests passing. Then, we should optimize. However, if you see
   something that could probably be optimized, let's not forget it. Add a // TODO: or // IMPROVEMENT: comment.

 - Coding conventions in this file (these may need to be changed before merging back into master)
   - Very basic hungarian notation:
     p prefix for pointers (pSelector)
     no prefix for value types and references (selector)
   - Use STL iterators where possible
   - prefer verbose naming over terse naming
   - use typedefs for STL container types for make maintenance easier

 - You may see a lot of comments that say "// TODO: is this the correct combinator?". See the comment referring to combinators
   in extendCompoundSelector for a more extensive explanation of my confusion. I think our divergence in data model from ruby
   sass causes this to be necessary.


 GLOBAL TODOS:

 - wrap the contents of the print functions in DEBUG preprocesser conditionals so they will be optimized away in non-debug mode.

 - consider making the extend* functions member functions to avoid passing around ctx and subset_map map around. This has the
   drawback that the implementation details of the operator are then exposed to the outside world, which is not ideal and
   can cause additional compile time dependencies.

 - mark the helper methods in this file static to given them compilation unit linkage.

 - implement parent directive matching

 - fix compilation warnings for unused Extend members if we really don't need those references anymore.
 */


namespace Sass {



#ifdef DEBUG

  // TODO: move the ast specific ostream operators into ast.hpp/ast.cpp
  std::ostream& operator<<(std::ostream& os, const Complex_Selector::Combinator combinator) {
    switch (combinator) {
      case Complex_Selector::ANCESTOR_OF: os << "\" \""; break;
      case Complex_Selector::PARENT_OF:   os << "\">\""; break;
      case Complex_Selector::PRECEDES:    os << "\"~\""; break;
      case Complex_Selector::ADJACENT_TO: os << "\"+\""; break;
      case Complex_Selector::REFERENCE:   os << "\"/\""; break;
    }

    return os;
  }


  std::ostream& operator<<(std::ostream& os, Compound_Selector& compoundSelector) {
    for (size_t i = 0, L = compoundSelector.length(); i < L; ++i) {
      if (i > 0) os << ", ";
      os << compoundSelector[i]->to_string();
    }
    return os;
  }

  std::ostream& operator<<(std::ostream& os, Simple_Selector& simpleSelector) {
    os << simpleSelector.to_string();
    return os;
  }

  // Print a string representation of a Compound_Selector
  static void printSimpleSelector(Simple_Selector* pSimpleSelector, const char* message=NULL, bool newline=true) {

    if (message) {
      std::cerr << message;
    }

    if (pSimpleSelector) {
      std::cerr << "[" << *pSimpleSelector << "]";
    } else {
      std::cerr << "NULL";
    }

    if (newline) {
      std::cerr << std::endl;
    }
  }

  // Print a string representation of a Compound_Selector
  static void printCompoundSelector(Compound_Selector_Ptr pCompoundSelector, const char* message=NULL, bool newline=true) {

    if (message) {
      std::cerr << message;
    }

    if (pCompoundSelector) {
      std::cerr << "[" << *pCompoundSelector << "]";
    } else {
      std::cerr << "NULL";
    }

    if (newline) {
      std::cerr << std::endl;
    }
  }


  std::ostream& operator<<(std::ostream& os, Complex_Selector& complexSelector) {

    os << "[";
    Complex_Selector_Ptr pIter = &complexSelector;
    bool first = true;
    while (pIter) {
      if (pIter->combinator() != Complex_Selector::ANCESTOR_OF) {
        if (!first) {
          os << ", ";
        }
        first = false;
        os << pIter->combinator();
      }

      if (!first) {
        os << ", ";
      }
      first = false;

      if (pIter->head()) {
        os << pIter->head()->to_string();
      } else {
        os << "NULL_HEAD";
      }

      pIter = pIter->tail();
    }
    os << "]";

    return os;
  }


  // Print a string representation of a Complex_Selector
  static void printComplexSelector(Complex_Selector_Ptr pComplexSelector, const char* message=NULL, bool newline=true) {

    if (message) {
      std::cerr << message;
    }

    if (pComplexSelector) {
      std::cerr << *pComplexSelector;
    } else {
      std::cerr << "NULL";
    }

    if (newline) {
      std::cerr << std::endl;
    }
  }

  static void printSelsNewSeqPairCollection(SubSetMapLookups& collection, const char* message=NULL, bool newline=true) {

    if (message) {
      std::cerr << message;
    }
    bool first = true;
    std::cerr << "[";
    for(SubSetMapLookup& pair : collection) {
      if (first) {
        first = false;
      } else {
        std::cerr << ", ";
      }
      std::cerr << "[";
      Compound_Selector_Ptr pSels = pair.first;
      Complex_Selector_Ptr pNewSelector = pair.second;
      std::cerr << "[" << *pSels << "], ";
      printComplexSelector(pNewSelector, NULL, false);
    }
    std::cerr << "]";

    if (newline) {
      std::cerr << std::endl;
    }
  }

  // Print a string representation of a ComplexSelectorSet
  static void printSourcesSet(ComplexSelectorSet& sources, const char* message=NULL, bool newline=true) {

    if (message) {
      std::cerr << message;
    }

    // Convert to a deque of strings so we can sort since order doesn't matter in a set. This should cut down on
    // the differences we see when debug printing.
    typedef std::deque<std::string> SourceStrings;
    SourceStrings sourceStrings;
    for (ComplexSelectorSet::iterator iterator = sources.begin(), iteratorEnd = sources.end(); iterator != iteratorEnd; ++iterator) {
      Complex_Selector_Ptr pSource = *iterator;
      std::stringstream sstream;
      sstream << complexSelectorToNode(pSource);
      sourceStrings.push_back(sstream.str());
    }

    // Sort to get consistent output
    std::sort(sourceStrings.begin(), sourceStrings.end());

    std::cerr << "ComplexSelectorSet[";
    for (SourceStrings::iterator iterator = sourceStrings.begin(), iteratorEnd = sourceStrings.end(); iterator != iteratorEnd; ++iterator) {
      std::string source = *iterator;
      if (iterator != sourceStrings.begin()) {
        std::cerr << ", ";
      }
      std::cerr << source;
    }
    std::cerr << "]";

    if (newline) {
      std::cerr << std::endl;
    }
  }


  std::ostream& operator<<(std::ostream& os, SubSetMapPairs& entries) {
    os << "SUBSET_MAP_ENTRIES[";

    for (SubSetMapPairs::iterator iterator = entries.begin(), endIterator = entries.end(); iterator != endIterator; ++iterator) {
      Complex_Selector_Obj pExtComplexSelector = iterator->first;    // The selector up to where the @extend is (ie, the thing to merge)
      Compound_Selector_Obj pExtCompoundSelector = iterator->second; // The stuff after the @extend

      if (iterator != entries.begin()) {
        os << ", ";
      }

      os << "(";

      if (pExtComplexSelector) {
        std::cerr << *pExtComplexSelector;
      } else {
        std::cerr << "NULL";
      }

      os << " -> ";

      if (pExtCompoundSelector) {
        std::cerr << *pExtCompoundSelector;
      } else {
        std::cerr << "NULL";
      }

      os << ")";

    }

    os << "]";

    return os;
  }
#endif

  static bool parentSuperselector(Complex_Selector_Ptr pOne, Complex_Selector_Ptr pTwo) {
    // TODO: figure out a better way to create a Complex_Selector from scratch
    // TODO: There's got to be a better way. This got ugly quick...
    Type_Selector_Obj fakeParent = SASS_MEMORY_NEW(Type_Selector, ParserState("[FAKE]"), "temp");
    Compound_Selector_Obj fakeHead = SASS_MEMORY_NEW(Compound_Selector, ParserState("[FAKE]"), 1 /*size*/);
    fakeHead->elements().push_back(fakeParent);
    Complex_Selector_Obj fakeParentContainer = SASS_MEMORY_NEW(Complex_Selector, ParserState("[FAKE]"), Complex_Selector::ANCESTOR_OF, fakeHead /*head*/, {} /*tail*/);

    pOne->set_innermost(fakeParentContainer, Complex_Selector::ANCESTOR_OF);
    pTwo->set_innermost(fakeParentContainer, Complex_Selector::ANCESTOR_OF);

    bool isSuperselector = pOne->is_superselector_of(pTwo);

    pOne->clear_innermost();
    pTwo->clear_innermost();

    return isSuperselector;
  }

  void nodeToComplexSelectorDeque(const Node& node, ComplexSelectorDeque& out) {
    for (NodeDeque::iterator iter = node.collection()->begin(), iterEnd = node.collection()->end(); iter != iterEnd; iter++) {
      Node& child = *iter;
      out.push_back(nodeToComplexSelector(child));
    }
  }

  Node complexSelectorDequeToNode(const ComplexSelectorDeque& deque) {
    Node result = Node::createCollection();

    for (ComplexSelectorDeque::const_iterator iter = deque.begin(), iterEnd = deque.end(); iter != iterEnd; iter++) {
      Complex_Selector_Obj pChild = *iter;
      result.collection()->push_back(complexSelectorToNode(pChild));
    }

    return result;
  }

  class LcsCollectionComparator {
  public:
    LcsCollectionComparator() {}

    bool operator()(Complex_Selector_Obj pOne, Complex_Selector_Obj pTwo, Complex_Selector_Obj& pOut) const {
      /*
      This code is based on the following block from ruby sass' subweave
        do |s1, s2|
          next s1 if s1 == s2
          next unless s1.first.is_a?(SimpleSequence) && s2.first.is_a?(SimpleSequence)
          next s2 if parent_superselector?(s1, s2)
          next s1 if parent_superselector?(s2, s1)
        end
      */

      if (*pOne == *pTwo) {
        pOut = pOne;
        return true;
      }

      if (pOne->combinator() != Complex_Selector::ANCESTOR_OF || pTwo->combinator() != Complex_Selector::ANCESTOR_OF) {
        return false;
      }

      if (parentSuperselector(pOne, pTwo)) {
        pOut = pTwo;
        return true;
      }

      if (parentSuperselector(pTwo, pOne)) {
        pOut = pOne;
        return true;
      }

      return false;
    }
  };


  /*
  This is the equivalent of ruby's Sass::Util.lcs_backtrace.

  # Computes a single longest common subsequence for arrays x and y.
  # Algorithm from http://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Reading_out_an_LCS
  */
  void lcs_backtrace(const LCSTable& c, ComplexSelectorDeque& x, ComplexSelectorDeque& y, int i, int j, const LcsCollectionComparator& comparator, ComplexSelectorDeque& out) {
    //DEBUG_PRINTLN(LCS, "LCSBACK: X=" << x << " Y=" << y << " I=" << i << " J=" << j)
    // TODO: make printComplexSelectorDeque and use DEBUG_EXEC AND DEBUG_PRINTLN HERE to get equivalent output

    if (i == 0 || j == 0) {
      DEBUG_PRINTLN(LCS, "RETURNING EMPTY")
      return;
    }


    Complex_Selector_Obj pCompareOut;
    if (comparator(x[i], y[j], pCompareOut)) {
      DEBUG_PRINTLN(LCS, "RETURNING AFTER ELEM COMPARE")
      lcs_backtrace(c, x, y, i - 1, j - 1, comparator, out);
      out.push_back(pCompareOut);
      return;
    }

    if (c[i][j - 1] > c[i - 1][j]) {
      DEBUG_PRINTLN(LCS, "RETURNING AFTER TABLE COMPARE")
      lcs_backtrace(c, x, y, i, j - 1, comparator, out);
      return;
    }

    DEBUG_PRINTLN(LCS, "FINAL RETURN")
    lcs_backtrace(c, x, y, i - 1, j, comparator, out);
    return;
  }

  /*
  This is the equivalent of ruby's Sass::Util.lcs_table.

  # Calculates the memoization table for the Least Common Subsequence algorithm.
  # Algorithm from http://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Computing_the_length_of_the_LCS
  */
  void lcs_table(const ComplexSelectorDeque& x, const ComplexSelectorDeque& y, const LcsCollectionComparator& comparator, LCSTable& out) {
    //DEBUG_PRINTLN(LCS, "LCSTABLE: X=" << x << " Y=" << y)
    // TODO: make printComplexSelectorDeque and use DEBUG_EXEC AND DEBUG_PRINTLN HERE to get equivalent output

    LCSTable c(x.size(), std::vector<int>(y.size()));

    // These shouldn't be necessary since the vector will be initialized to 0 already.
    // x.size.times {|i| c[i][0] = 0}
    // y.size.times {|j| c[0][j] = 0}

    for (size_t i = 1; i < x.size(); i++) {
      for (size_t j = 1; j < y.size(); j++) {
        Complex_Selector_Obj pCompareOut;

        if (comparator(x[i], y[j], pCompareOut)) {
          c[i][j] = c[i - 1][j - 1] + 1;
        } else {
          c[i][j] = std::max(c[i][j - 1], c[i - 1][j]);
        }
      }
    }

    out = c;
  }

  /*
  This is the equivalent of ruby's Sass::Util.lcs.

  # Computes a single longest common subsequence for `x` and `y`.
  # If there are more than one longest common subsequences,
  # the one returned is that which starts first in `x`.

  # @param x [NodeCollection]
  # @param y [NodeCollection]
  # @comparator An equality check between elements of `x` and `y`.
  # @return [NodeCollection] The LCS

  http://en.wikipedia.org/wiki/Longest_common_subsequence_problem
  */
  void lcs(ComplexSelectorDeque& x, ComplexSelectorDeque& y, const LcsCollectionComparator& comparator, ComplexSelectorDeque& out) {
    //DEBUG_PRINTLN(LCS, "LCS: X=" << x << " Y=" << y)
    // TODO: make printComplexSelectorDeque and use DEBUG_EXEC AND DEBUG_PRINTLN HERE to get equivalent output

    x.push_front({});
    y.push_front({});

    LCSTable table;
    lcs_table(x, y, comparator, table);

    return lcs_backtrace(table, x, y, static_cast<int>(x.size()) - 1, static_cast<int>(y.size()) - 1, comparator, out);
  }


  /*
   This is the equivalent of ruby's Sequence.trim.

   The following is the modified version of the ruby code that was more portable to C++. You
   should be able to drop it into ruby 3.2.19 and get the same results from ruby sass.

        # Avoid truly horrific quadratic behavior. TODO: I think there
        # may be a way to get perfect trimming without going quadratic.
        return seqses if seqses.size > 100

        # Keep the results in a separate array so we can be sure we aren't
        # comparing against an already-trimmed selector. This ensures that two
        # identical selectors don't mutually trim one another.
        result = seqses.dup

        # This is n^2 on the sequences, but only comparing between
        # separate sequences should limit the quadratic behavior.
        seqses.each_with_index do |seqs1, i|
          tempResult = []

          for seq1 in seqs1 do
            max_spec = 0
            for seq in _sources(seq1) do
              max_spec = [max_spec, seq.specificity].max
            end


            isMoreSpecificOuter = false
            for seqs2 in result do
              if seqs1.equal?(seqs2) then
                next
              end

              # Second Law of Extend: the specificity of a generated selector
              # should never be less than the specificity of the extending
              # selector.
              #
              # See https://github.com/nex3/sass/issues/324.
              isMoreSpecificInner = false
              for seq2 in seqs2 do
                isMoreSpecificInner = _specificity(seq2) >= max_spec && _superselector?(seq2, seq1)
                if isMoreSpecificInner then
                  break
                end
              end

              if isMoreSpecificInner then
                isMoreSpecificOuter = true
                break
              end
            end

            if !isMoreSpecificOuter then
              tempResult.push(seq1)
            end
          end

          result[i] = tempResult

        end

        result
   */
  /*
   - IMPROVEMENT: We could probably work directly in the output trimmed deque.
   */
  Node Extend::trim(Node& seqses, bool isReplace) {
    // See the comments in the above ruby code before embarking on understanding this function.

    // Avoid poor performance in extreme cases.
    if (seqses.collection()->size() > 100) {
      return seqses;
    }


    DEBUG_PRINTLN(TRIM, "TRIM: " << seqses)


    Node result = Node::createCollection();
    result.plus(seqses);

    DEBUG_PRINTLN(TRIM, "RESULT INITIAL: " << result)

    // Normally we use the standard STL iterators, but in this case, we need to access the result collection by index since we're
    // iterating the input collection, computing a value, and then setting the result in the output collection. We have to keep track
    // of the index manually.
    int toTrimIndex = 0;

    for (NodeDeque::iterator seqsesIter = seqses.collection()->begin(), seqsesIterEnd = seqses.collection()->end(); seqsesIter != seqsesIterEnd; ++seqsesIter) {
      Node& seqs1 = *seqsesIter;

      DEBUG_PRINTLN(TRIM, "SEQS1: " << seqs1 << " " << toTrimIndex)

      Node tempResult = Node::createCollection();
      tempResult.got_line_feed = seqs1.got_line_feed;

      for (NodeDeque::iterator seqs1Iter = seqs1.collection()->begin(), seqs1EndIter = seqs1.collection()->end(); seqs1Iter != seqs1EndIter; ++seqs1Iter) {
        Node& seq1 = *seqs1Iter;

        Complex_Selector_Obj pSeq1 = nodeToComplexSelector(seq1);

        // Compute the maximum specificity. This requires looking at the "sources" of the sequence. See SimpleSequence.sources in the ruby code
        // for a good description of sources.
        //
        // TODO: I'm pretty sure there's a bug in the sources code. It was implemented for sass-spec's 182_test_nested_extend_loop test.
        // While the test passes, I compared the state of each trim call to verify correctness. The last trim call had incorrect sources. We
        // had an extra source that the ruby version did not have. Without a failing test case, this is going to be extra hard to find. My
        // best guess at this point is that we're cloning an object somewhere and maintaining the sources when we shouldn't be. This is purely
        // a guess though.
        unsigned long maxSpecificity = isReplace ? pSeq1->specificity() : 0;
        ComplexSelectorSet sources = pSeq1->sources();

        DEBUG_PRINTLN(TRIM, "TRIM SEQ1: " << seq1)
        DEBUG_EXEC(TRIM, printSourcesSet(sources, "TRIM SOURCES: "))

        for (ComplexSelectorSet::iterator sourcesSetIterator = sources.begin(), sourcesSetIteratorEnd = sources.end(); sourcesSetIterator != sourcesSetIteratorEnd; ++sourcesSetIterator) {
          const Complex_Selector_Obj& pCurrentSelector = *sourcesSetIterator;
          maxSpecificity = std::max(maxSpecificity, pCurrentSelector->specificity());
        }

        DEBUG_PRINTLN(TRIM, "MAX SPECIFICITY: " << maxSpecificity)

        bool isMoreSpecificOuter = false;

        int resultIndex = 0;

        for (NodeDeque::iterator resultIter = result.collection()->begin(), resultIterEnd = result.collection()->end(); resultIter != resultIterEnd; ++resultIter) {
          Node& seqs2 = *resultIter;

          DEBUG_PRINTLN(TRIM, "SEQS1: " << seqs1)
          DEBUG_PRINTLN(TRIM, "SEQS2: " << seqs2)

          // Do not compare the same sequence to itself. The ruby call we're trying to
          // emulate is: seqs1.equal?(seqs2). equal? is an object comparison, not an equivalency comparision.
          // Since we have the same pointers in seqes and results, we can do a pointer comparision. seqs1 is
          // derived from seqses and seqs2 is derived from result.
          if (seqs1.collection() == seqs2.collection()) {
            DEBUG_PRINTLN(TRIM, "CONTINUE")
            continue;
          }

          bool isMoreSpecificInner = false;

          for (NodeDeque::iterator seqs2Iter = seqs2.collection()->begin(), seqs2IterEnd = seqs2.collection()->end(); seqs2Iter != seqs2IterEnd; ++seqs2Iter) {
            Node& seq2 = *seqs2Iter;

            Complex_Selector_Obj pSeq2 = nodeToComplexSelector(seq2);

            DEBUG_PRINTLN(TRIM, "SEQ2 SPEC: " << pSeq2->specificity())
            DEBUG_PRINTLN(TRIM, "IS SPEC: " << pSeq2->specificity() << " >= " << maxSpecificity << " " << (pSeq2->specificity() >= maxSpecificity ? "true" : "false"))
            DEBUG_PRINTLN(TRIM, "IS SUPER: " << (pSeq2->is_superselector_of(pSeq1) ? "true" : "false"))

            isMoreSpecificInner = pSeq2->specificity() >= maxSpecificity && pSeq2->is_superselector_of(pSeq1);

            if (isMoreSpecificInner) {
              DEBUG_PRINTLN(TRIM, "FOUND MORE SPECIFIC")
              break;
            }
          }

          // If we found something more specific, we're done. Let the outer loop know and stop iterating.
          if (isMoreSpecificInner) {
            isMoreSpecificOuter = true;
            break;
          }

          resultIndex++;
        }

        if (!isMoreSpecificOuter) {
          DEBUG_PRINTLN(TRIM, "PUSHING: " << seq1)
          tempResult.collection()->push_back(seq1);
        }

      }

      DEBUG_PRINTLN(TRIM, "RESULT BEFORE ASSIGN: " << result)
      DEBUG_PRINTLN(TRIM, "TEMP RESULT: " << toTrimIndex << " " << tempResult)
      (*result.collection())[toTrimIndex] = tempResult;

      toTrimIndex++;

      DEBUG_PRINTLN(TRIM, "RESULT: " << result)
    }

    return result;
  }



  static bool parentSuperselector(const Node& one, const Node& two) {
    // TODO: figure out a better way to create a Complex_Selector from scratch
    // TODO: There's got to be a better way. This got ugly quick...
    Type_Selector_Obj fakeParent = SASS_MEMORY_NEW(Type_Selector, ParserState("[FAKE]"), "temp");
    Compound_Selector_Obj fakeHead = SASS_MEMORY_NEW(Compound_Selector, ParserState("[FAKE]"), 1 /*size*/);
    fakeHead->elements().push_back(fakeParent);
    Complex_Selector_Obj fakeParentContainer = SASS_MEMORY_NEW(Complex_Selector, ParserState("[FAKE]"), Complex_Selector::ANCESTOR_OF, fakeHead /*head*/, {} /*tail*/);

    Complex_Selector_Obj pOneWithFakeParent = nodeToComplexSelector(one);
    pOneWithFakeParent->set_innermost(fakeParentContainer, Complex_Selector::ANCESTOR_OF);
    Complex_Selector_Obj pTwoWithFakeParent = nodeToComplexSelector(two);
    pTwoWithFakeParent->set_innermost(fakeParentContainer, Complex_Selector::ANCESTOR_OF);

    return pOneWithFakeParent->is_superselector_of(pTwoWithFakeParent);
  }


  class ParentSuperselectorChunker {
  public:
    ParentSuperselectorChunker(Node& lcs) : mLcs(lcs) {}
    Node& mLcs;

    bool operator()(const Node& seq) const {
      // {|s| parent_superselector?(s.first, lcs.first)}
      if (seq.collection()->size() == 0) return false;
      return parentSuperselector(seq.collection()->front(), mLcs.collection()->front());
    }
  };

  class SubweaveEmptyChunker {
  public:
    bool operator()(const Node& seq) const {
      // {|s| s.empty?}

      return seq.collection()->empty();
    }
  };

  /*
  # Takes initial subsequences of `seq1` and `seq2` and returns all
  # orderings of those subsequences. The initial subsequences are determined
  # by a block.
  #
  # Destructively removes the initial subsequences of `seq1` and `seq2`.
  #
  # For example, given `(A B C | D E)` and `(1 2 | 3 4 5)` (with `|`
  # denoting the boundary of the initial subsequence), this would return
  # `[(A B C 1 2), (1 2 A B C)]`. The sequences would then be `(D E)` and
  # `(3 4 5)`.
  #
  # @param seq1 [Array]
  # @param seq2 [Array]
  # @yield [a] Used to determine when to cut off the initial subsequences.
  #   Called repeatedly for each sequence until it returns true.
  # @yieldparam a [Array] A final subsequence of one input sequence after
  #   cutting off some initial subsequence.
  # @yieldreturn [Boolean] Whether or not to cut off the initial subsequence
  #   here.
  # @return [Array<Array>] All possible orderings of the initial subsequences.
  def chunks(seq1, seq2)
    chunk1 = []
    chunk1 << seq1.shift until yield seq1
    chunk2 = []
    chunk2 << seq2.shift until yield seq2
    return [] if chunk1.empty? && chunk2.empty?
    return [chunk2] if chunk1.empty?
    return [chunk1] if chunk2.empty?
    [chunk1 + chunk2, chunk2 + chunk1]
  end
  */
  template<typename ChunkerType>
  static Node chunks(Node& seq1, Node& seq2, const ChunkerType& chunker) {
    Node chunk1 = Node::createCollection();
    while (seq1.collection()->size() && !chunker(seq1)) {
      chunk1.collection()->push_back(seq1.collection()->front());
      seq1.collection()->pop_front();
    }

    Node chunk2 = Node::createCollection();
    while (!seq2.collection()->empty() && !chunker(seq2)) {
      chunk2.collection()->push_back(seq2.collection()->front());
      seq2.collection()->pop_front();
    }

    if (chunk1.collection()->empty() && chunk2.collection()->empty()) {
      DEBUG_PRINTLN(CHUNKS, "RETURNING BOTH EMPTY")
      return Node::createCollection();
    }

    if (chunk1.collection()->empty()) {
      Node chunk2Wrapper = Node::createCollection();
      chunk2Wrapper.collection()->push_back(chunk2);
      DEBUG_PRINTLN(CHUNKS, "RETURNING ONE EMPTY")
      return chunk2Wrapper;
    }

    if (chunk2.collection()->empty()) {
      Node chunk1Wrapper = Node::createCollection();
      chunk1Wrapper.collection()->push_back(chunk1);
      DEBUG_PRINTLN(CHUNKS, "RETURNING TWO EMPTY")
      return chunk1Wrapper;
    }

    Node perms = Node::createCollection();

    Node firstPermutation = Node::createCollection();
    firstPermutation.collection()->insert(firstPermutation.collection()->end(), chunk1.collection()->begin(), chunk1.collection()->end());
    firstPermutation.collection()->insert(firstPermutation.collection()->end(), chunk2.collection()->begin(), chunk2.collection()->end());
    perms.collection()->push_back(firstPermutation);

    Node secondPermutation = Node::createCollection();
    secondPermutation.collection()->insert(secondPermutation.collection()->end(), chunk2.collection()->begin(), chunk2.collection()->end());
    secondPermutation.collection()->insert(secondPermutation.collection()->end(), chunk1.collection()->begin(), chunk1.collection()->end());
    perms.collection()->push_back(secondPermutation);

    DEBUG_PRINTLN(CHUNKS, "RETURNING PERM")

    return perms;
  }


  static Node groupSelectors(Node& seq) {
    Node newSeq = Node::createCollection();

    Node tail = Node::createCollection();
    tail.plus(seq);

    while (!tail.collection()->empty()) {
      Node head = Node::createCollection();

      do {
        head.collection()->push_back(tail.collection()->front());
        tail.collection()->pop_front();
      } while (!tail.collection()->empty() && (head.collection()->back().isCombinator() || tail.collection()->front().isCombinator()));

      newSeq.collection()->push_back(head);
    }

    return newSeq;
  }


  static void getAndRemoveInitialOps(Node& seq, Node& ops) {
    NodeDeque& seqCollection = *(seq.collection());
    NodeDeque& opsCollection = *(ops.collection());

    while (seqCollection.size() > 0 && seqCollection.front().isCombinator()) {
      opsCollection.push_back(seqCollection.front());
      seqCollection.pop_front();
    }
  }


  static void getAndRemoveFinalOps(Node& seq, Node& ops) {
    NodeDeque& seqCollection = *(seq.collection());
    NodeDeque& opsCollection = *(ops.collection());

    while (seqCollection.size() > 0 && seqCollection.back().isCombinator()) {
      opsCollection.push_back(seqCollection.back()); // Purposefully reversed to match ruby code
      seqCollection.pop_back();
    }
  }


  /*
      def merge_initial_ops(seq1, seq2)
        ops1, ops2 = [], []
        ops1 << seq1.shift while seq1.first.is_a?(String)
        ops2 << seq2.shift while seq2.first.is_a?(String)

        newline = false
        newline ||= !!ops1.shift if ops1.first == "\n"
        newline ||= !!ops2.shift if ops2.first == "\n"

        # If neither sequence is a subsequence of the other, they cannot be
        # merged successfully
        lcs = Sass::Util.lcs(ops1, ops2)
        return unless lcs == ops1 || lcs == ops2
        return (newline ? ["\n"] : []) + (ops1.size > ops2.size ? ops1 : ops2)
      end
  */
  static Node mergeInitialOps(Node& seq1, Node& seq2) {
    Node ops1 = Node::createCollection();
    Node ops2 = Node::createCollection();

    getAndRemoveInitialOps(seq1, ops1);
    getAndRemoveInitialOps(seq2, ops2);

    // TODO: Do we have this information available to us?
    // newline = false
    // newline ||= !!ops1.shift if ops1.first == "\n"
    // newline ||= !!ops2.shift if ops2.first == "\n"

    // If neither sequence is a subsequence of the other, they cannot be merged successfully
    DefaultLcsComparator lcsDefaultComparator;
    Node opsLcs = lcs(ops1, ops2, lcsDefaultComparator);

    if (!(opsLcs == ops1 || opsLcs == ops2)) {
      return Node::createNil();
    }

    // TODO: more newline logic
    // return (newline ? ["\n"] : []) + (ops1.size > ops2.size ? ops1 : ops2)

    return (ops1.collection()->size() > ops2.collection()->size() ? ops1 : ops2);
  }


  /*
      def merge_final_ops(seq1, seq2, res = [])


        # This code looks complicated, but it's actually just a bunch of special
        # cases for interactions between different combinators.
        op1, op2 = ops1.first, ops2.first
        if op1 && op2
          sel1 = seq1.pop
          sel2 = seq2.pop
          if op1 == '~' && op2 == '~'
            if sel1.superselector?(sel2)
              res.unshift sel2, '~'
            elsif sel2.superselector?(sel1)
              res.unshift sel1, '~'
            else
              merged = sel1.unify(sel2.members, sel2.subject?)
              res.unshift [
                [sel1, '~', sel2, '~'],
                [sel2, '~', sel1, '~'],
                ([merged, '~'] if merged)
              ].compact
            end
          elsif (op1 == '~' && op2 == '+') || (op1 == '+' && op2 == '~')
            if op1 == '~'
              tilde_sel, plus_sel = sel1, sel2
            else
              tilde_sel, plus_sel = sel2, sel1
            end

            if tilde_sel.superselector?(plus_sel)
              res.unshift plus_sel, '+'
            else
              merged = plus_sel.unify(tilde_sel.members, tilde_sel.subject?)
              res.unshift [
                [tilde_sel, '~', plus_sel, '+'],
                ([merged, '+'] if merged)
              ].compact
            end
          elsif op1 == '>' && %w[~ +].include?(op2)
            res.unshift sel2, op2
            seq1.push sel1, op1
          elsif op2 == '>' && %w[~ +].include?(op1)
            res.unshift sel1, op1
            seq2.push sel2, op2
          elsif op1 == op2
            return unless merged = sel1.unify(sel2.members, sel2.subject?)
            res.unshift merged, op1
          else
            # Unknown selector combinators can't be unified
            return
          end
          return merge_final_ops(seq1, seq2, res)
        elsif op1
          seq2.pop if op1 == '>' && seq2.last && seq2.last.superselector?(seq1.last)
          res.unshift seq1.pop, op1
          return merge_final_ops(seq1, seq2, res)
        else # op2
          seq1.pop if op2 == '>' && seq1.last && seq1.last.superselector?(seq2.last)
          res.unshift seq2.pop, op2
          return merge_final_ops(seq1, seq2, res)
        end
      end
  */
  static Node mergeFinalOps(Node& seq1, Node& seq2, Node& res) {

    Node ops1 = Node::createCollection();
    Node ops2 = Node::createCollection();

    getAndRemoveFinalOps(seq1, ops1);
    getAndRemoveFinalOps(seq2, ops2);

    // TODO: do we have newlines to remove?
    // ops1.reject! {|o| o == "\n"}
    // ops2.reject! {|o| o == "\n"}

    if (ops1.collection()->empty() && ops2.collection()->empty()) {
      return res;
    }

    if (ops1.collection()->size() > 1 || ops2.collection()->size() > 1) {
      DefaultLcsComparator lcsDefaultComparator;
      Node opsLcs = lcs(ops1, ops2, lcsDefaultComparator);

      // If there are multiple operators, something hacky's going on. If one is a supersequence of the other, use that, otherwise give up.

      if (!(opsLcs == ops1 || opsLcs == ops2)) {
        return Node::createNil();
      }

      if (ops1.collection()->size() > ops2.collection()->size()) {
        res.collection()->insert(res.collection()->begin(), ops1.collection()->rbegin(), ops1.collection()->rend());
      } else {
        res.collection()->insert(res.collection()->begin(), ops2.collection()->rbegin(), ops2.collection()->rend());
      }

      return res;
    }

    if (!ops1.collection()->empty() && !ops2.collection()->empty()) {

      Node op1 = ops1.collection()->front();
      Node op2 = ops2.collection()->front();

      Node sel1 = seq1.collection()->back();
      seq1.collection()->pop_back();

      Node sel2 = seq2.collection()->back();
      seq2.collection()->pop_back();

      if (op1.combinator() == Complex_Selector::PRECEDES && op2.combinator() == Complex_Selector::PRECEDES) {

        if (sel1.selector()->is_superselector_of(sel2.selector())) {

          res.collection()->push_front(op1 /*PRECEDES - could have been op2 as well*/);
          res.collection()->push_front(sel2);

        } else if (sel2.selector()->is_superselector_of(sel1.selector())) {

          res.collection()->push_front(op1 /*PRECEDES - could have been op2 as well*/);
          res.collection()->push_front(sel1);

        } else {

          DEBUG_PRINTLN(ALL, "sel1: " << sel1)
          DEBUG_PRINTLN(ALL, "sel2: " << sel2)

          Complex_Selector_Obj pMergedWrapper = SASS_MEMORY_CLONE(sel1.selector()); // Clone the Complex_Selector to get back to something we can transform to a node once we replace the head with the unification result
          // TODO: does subject matter? Ruby: return unless merged = sel1.unify(sel2.members, sel2.subject?)
          Compound_Selector_Ptr pMerged = sel1.selector()->head()->unify_with(sel2.selector()->head());
          pMergedWrapper->head(pMerged);

          DEBUG_EXEC(ALL, printCompoundSelector(pMerged, "MERGED: "))

          Node newRes = Node::createCollection();

          Node firstPerm = Node::createCollection();
          firstPerm.collection()->push_back(sel1);
          firstPerm.collection()->push_back(Node::createCombinator(Complex_Selector::PRECEDES));
          firstPerm.collection()->push_back(sel2);
          firstPerm.collection()->push_back(Node::createCombinator(Complex_Selector::PRECEDES));
          newRes.collection()->push_back(firstPerm);

          Node secondPerm = Node::createCollection();
          secondPerm.collection()->push_back(sel2);
          secondPerm.collection()->push_back(Node::createCombinator(Complex_Selector::PRECEDES));
          secondPerm.collection()->push_back(sel1);
          secondPerm.collection()->push_back(Node::createCombinator(Complex_Selector::PRECEDES));
          newRes.collection()->push_back(secondPerm);

          if (pMerged) {
            Node mergedPerm = Node::createCollection();
            mergedPerm.collection()->push_back(Node::createSelector(pMergedWrapper));
            mergedPerm.collection()->push_back(Node::createCombinator(Complex_Selector::PRECEDES));
            newRes.collection()->push_back(mergedPerm);
          }

          res.collection()->push_front(newRes);

          DEBUG_PRINTLN(ALL, "RESULT: " << res)

        }

      } else if (((op1.combinator() == Complex_Selector::PRECEDES && op2.combinator() == Complex_Selector::ADJACENT_TO)) || ((op1.combinator() == Complex_Selector::ADJACENT_TO && op2.combinator() == Complex_Selector::PRECEDES))) {

          Node tildeSel = sel1;
          Node plusSel = sel2;
          Node plusOp = op2;
          if (op1.combinator() != Complex_Selector::PRECEDES) {
            tildeSel = sel2;
            plusSel = sel1;
            plusOp = op1;
          }

          if (tildeSel.selector()->is_superselector_of(plusSel.selector())) {

            res.collection()->push_front(plusOp);
            res.collection()->push_front(plusSel);

          } else {

            DEBUG_PRINTLN(ALL, "PLUS SEL: " << plusSel)
            DEBUG_PRINTLN(ALL, "TILDE SEL: " << tildeSel)

            Complex_Selector_Obj pMergedWrapper = SASS_MEMORY_CLONE(plusSel.selector()); // Clone the Complex_Selector to get back to something we can transform to a node once we replace the head with the unification result
            // TODO: does subject matter? Ruby: merged = plus_sel.unify(tilde_sel.members, tilde_sel.subject?)
            Compound_Selector_Ptr pMerged = plusSel.selector()->head()->unify_with(tildeSel.selector()->head());
            pMergedWrapper->head(pMerged);

            DEBUG_EXEC(ALL, printCompoundSelector(pMerged, "MERGED: "))

            Node newRes = Node::createCollection();

            Node firstPerm = Node::createCollection();
            firstPerm.collection()->push_back(tildeSel);
            firstPerm.collection()->push_back(Node::createCombinator(Complex_Selector::PRECEDES));
            firstPerm.collection()->push_back(plusSel);
            firstPerm.collection()->push_back(Node::createCombinator(Complex_Selector::ADJACENT_TO));
            newRes.collection()->push_back(firstPerm);

            if (pMerged) {
              Node mergedPerm = Node::createCollection();
              mergedPerm.collection()->push_back(Node::createSelector(pMergedWrapper));
              mergedPerm.collection()->push_back(Node::createCombinator(Complex_Selector::ADJACENT_TO));
              newRes.collection()->push_back(mergedPerm);
            }

            res.collection()->push_front(newRes);

            DEBUG_PRINTLN(ALL, "RESULT: " << res)

          }
      } else if (op1.combinator() == Complex_Selector::PARENT_OF && (op2.combinator() == Complex_Selector::PRECEDES || op2.combinator() == Complex_Selector::ADJACENT_TO)) {

        res.collection()->push_front(op2);
        res.collection()->push_front(sel2);

        seq1.collection()->push_back(sel1);
        seq1.collection()->push_back(op1);

      } else if (op2.combinator() == Complex_Selector::PARENT_OF && (op1.combinator() == Complex_Selector::PRECEDES || op1.combinator() == Complex_Selector::ADJACENT_TO)) {

        res.collection()->push_front(op1);
        res.collection()->push_front(sel1);

        seq2.collection()->push_back(sel2);
        seq2.collection()->push_back(op2);

      } else if (op1.combinator() == op2.combinator()) {

        DEBUG_PRINTLN(ALL, "sel1: " << sel1)
        DEBUG_PRINTLN(ALL, "sel2: " << sel2)

        Complex_Selector_Obj pMergedWrapper = SASS_MEMORY_CLONE(sel1.selector()); // Clone the Complex_Selector to get back to something we can transform to a node once we replace the head with the unification result
        // TODO: does subject matter? Ruby: return unless merged = sel1.unify(sel2.members, sel2.subject?)
        Compound_Selector_Ptr pMerged = sel1.selector()->head()->unify_with(sel2.selector()->head());
        pMergedWrapper->head(pMerged);

        DEBUG_EXEC(ALL, printCompoundSelector(pMerged, "MERGED: "))

        if (!pMerged) {
          return Node::createNil();
        }

        res.collection()->push_front(op1);
        res.collection()->push_front(Node::createSelector(pMergedWrapper));

        DEBUG_PRINTLN(ALL, "RESULT: " << res)

      } else {
        return Node::createNil();
      }

      return mergeFinalOps(seq1, seq2, res);

    } else if (!ops1.collection()->empty()) {

      Node op1 = ops1.collection()->front();

      if (op1.combinator() == Complex_Selector::PARENT_OF && !seq2.collection()->empty() && seq2.collection()->back().selector()->is_superselector_of(seq1.collection()->back().selector())) {
        seq2.collection()->pop_back();
      }

      // TODO: consider unshift(NodeCollection, Node)
      res.collection()->push_front(op1);
      res.collection()->push_front(seq1.collection()->back());
      seq1.collection()->pop_back();

      return mergeFinalOps(seq1, seq2, res);

    } else { // !ops2.collection()->empty()

      Node op2 = ops2.collection()->front();

      if (op2.combinator() == Complex_Selector::PARENT_OF && !seq1.collection()->empty() && seq1.collection()->back().selector()->is_superselector_of(seq2.collection()->back().selector())) {
        seq1.collection()->pop_back();
      }

      res.collection()->push_front(op2);
      res.collection()->push_front(seq2.collection()->back());
      seq2.collection()->pop_back();

      return mergeFinalOps(seq1, seq2, res);

    }

  }


  /*
    This is the equivalent of ruby's Sequence.subweave.

    Here is the original subweave code for reference during porting.

      def subweave(seq1, seq2)
        return [seq2] if seq1.empty?
        return [seq1] if seq2.empty?

        seq1, seq2 = seq1.dup, seq2.dup
        return unless init = merge_initial_ops(seq1, seq2)
        return unless fin = merge_final_ops(seq1, seq2)
        seq1 = group_selectors(seq1)
        seq2 = group_selectors(seq2)
        lcs = Sass::Util.lcs(seq2, seq1) do |s1, s2|
          next s1 if s1 == s2
          next unless s1.first.is_a?(SimpleSequence) && s2.first.is_a?(SimpleSequence)
          next s2 if parent_superselector?(s1, s2)
          next s1 if parent_superselector?(s2, s1)
        end

        diff = [[init]]
        until lcs.empty?
          diff << chunks(seq1, seq2) {|s| parent_superselector?(s.first, lcs.first)} << [lcs.shift]
          seq1.shift
          seq2.shift
        end
        diff << chunks(seq1, seq2) {|s| s.empty?}
        diff += fin.map {|sel| sel.is_a?(Array) ? sel : [sel]}
        diff.reject! {|c| c.empty?}

        result = Sass::Util.paths(diff).map {|p| p.flatten}.reject {|p| path_has_two_subjects?(p)}

        result
      end
  */
  Node subweave(Node& one, Node& two) {
    // Check for the simple cases
    if (one.collection()->size() == 0) {
      Node out = Node::createCollection();
      out.collection()->push_back(two);
      return out;
    }
    if (two.collection()->size() == 0) {
      Node out = Node::createCollection();
      out.collection()->push_back(one);
      return out;
    }

    Node seq1 = Node::createCollection();
    seq1.plus(one);
    Node seq2 = Node::createCollection();
    seq2.plus(two);

    DEBUG_PRINTLN(SUBWEAVE, "SUBWEAVE ONE: " << seq1)
    DEBUG_PRINTLN(SUBWEAVE, "SUBWEAVE TWO: " << seq2)

    Node init = mergeInitialOps(seq1, seq2);
    if (init.isNil()) {
      return Node::createNil();
    }

    DEBUG_PRINTLN(SUBWEAVE, "INIT: " << init)

    Node res = Node::createCollection();
    Node fin = mergeFinalOps(seq1, seq2, res);
    if (fin.isNil()) {
      return Node::createNil();
    }

    DEBUG_PRINTLN(SUBWEAVE, "FIN: " << fin)


    // Moving this line up since fin isn't modified between now and when it happened before
    // fin.map {|sel| sel.is_a?(Array) ? sel : [sel]}

    for (NodeDeque::iterator finIter = fin.collection()->begin(), finEndIter = fin.collection()->end();
           finIter != finEndIter; ++finIter) {

      Node& childNode = *finIter;

      if (!childNode.isCollection()) {
        Node wrapper = Node::createCollection();
        wrapper.collection()->push_back(childNode);
        childNode = wrapper;
      }

    }

    DEBUG_PRINTLN(SUBWEAVE, "FIN MAPPED: " << fin)



    Node groupSeq1 = groupSelectors(seq1);
    DEBUG_PRINTLN(SUBWEAVE, "SEQ1: " << groupSeq1)

    Node groupSeq2 = groupSelectors(seq2);
    DEBUG_PRINTLN(SUBWEAVE, "SEQ2: " << groupSeq2)


    ComplexSelectorDeque groupSeq1Converted;
    nodeToComplexSelectorDeque(groupSeq1, groupSeq1Converted);

    ComplexSelectorDeque groupSeq2Converted;
    nodeToComplexSelectorDeque(groupSeq2, groupSeq2Converted);

    ComplexSelectorDeque out;
    LcsCollectionComparator collectionComparator;
    lcs(groupSeq2Converted, groupSeq1Converted, collectionComparator, out);
    Node seqLcs = complexSelectorDequeToNode(out);

    DEBUG_PRINTLN(SUBWEAVE, "SEQLCS: " << seqLcs)


    Node initWrapper = Node::createCollection();
    initWrapper.collection()->push_back(init);
    Node diff = Node::createCollection();
    diff.collection()->push_back(initWrapper);

    DEBUG_PRINTLN(SUBWEAVE, "DIFF INIT: " << diff)


    while (!seqLcs.collection()->empty()) {
      ParentSuperselectorChunker superselectorChunker(seqLcs);
      Node chunksResult = chunks(groupSeq1, groupSeq2, superselectorChunker);
      diff.collection()->push_back(chunksResult);

      Node lcsWrapper = Node::createCollection();
      lcsWrapper.collection()->push_back(seqLcs.collection()->front());
      seqLcs.collection()->pop_front();
      diff.collection()->push_back(lcsWrapper);

      if (groupSeq1.collection()->size()) groupSeq1.collection()->pop_front();
      if (groupSeq2.collection()->size()) groupSeq2.collection()->pop_front();
    }

    DEBUG_PRINTLN(SUBWEAVE, "DIFF POST LCS: " << diff)


    DEBUG_PRINTLN(SUBWEAVE, "CHUNKS: ONE=" << groupSeq1 << " TWO=" << groupSeq2)


    SubweaveEmptyChunker emptyChunker;
    Node chunksResult = chunks(groupSeq1, groupSeq2, emptyChunker);
    diff.collection()->push_back(chunksResult);


    DEBUG_PRINTLN(SUBWEAVE, "DIFF POST CHUNKS: " << diff)


    diff.collection()->insert(diff.collection()->end(), fin.collection()->begin(), fin.collection()->end());

    DEBUG_PRINTLN(SUBWEAVE, "DIFF POST FIN MAPPED: " << diff)

    // JMA - filter out the empty nodes (use a new collection, since iterator erase() invalidates the old collection)
    Node diffFiltered = Node::createCollection();
    for (NodeDeque::iterator diffIter = diff.collection()->begin(), diffEndIter = diff.collection()->end();
           diffIter != diffEndIter; ++diffIter) {
      Node& node = *diffIter;
      if (node.collection() && !node.collection()->empty()) {
        diffFiltered.collection()->push_back(node);
      }
    }
    diff = diffFiltered;

    DEBUG_PRINTLN(SUBWEAVE, "DIFF POST REJECT: " << diff)


    Node pathsResult = paths(diff);

    DEBUG_PRINTLN(SUBWEAVE, "PATHS: " << pathsResult)


    // We're flattening in place
    for (NodeDeque::iterator pathsIter = pathsResult.collection()->begin(), pathsEndIter = pathsResult.collection()->end();
      pathsIter != pathsEndIter; ++pathsIter) {

      Node& child = *pathsIter;
      child = flatten(child);
    }

    DEBUG_PRINTLN(SUBWEAVE, "FLATTENED: " << pathsResult)


    /*
      TODO: implement
      rejected = mapped.reject {|p| path_has_two_subjects?(p)}
      $stderr.puts "REJECTED: #{rejected}"
     */


    return pathsResult;

  }
  /*
  // disabled to avoid clang warning [-Wunused-function]
  static Node subweaveNaive(const Node& one, const Node& two) {
    Node out = Node::createCollection();

    // Check for the simple cases
    if (one.isNil()) {
      out.collection()->push_back(two.klone());
    } else if (two.isNil()) {
      out.collection()->push_back(one.klone());
    } else {
      // Do the naive implementation. pOne = A B and pTwo = C D ...yields...  A B C D and C D A B
      // See https://gist.github.com/nex3/7609394 for details.

      Node firstPerm = one.klone();
      Node twoCloned = two.klone();
      firstPerm.plus(twoCloned);
      out.collection()->push_back(firstPerm);

      Node secondPerm = two.klone();
      Node oneCloned = one.klone();
      secondPerm.plus(oneCloned );
      out.collection()->push_back(secondPerm);
    }

    return out;
  }
  */


  /*
   This is the equivalent of ruby's Sequence.weave.

   The following is the modified version of the ruby code that was more portable to C++. You
   should be able to drop it into ruby 3.2.19 and get the same results from ruby sass.

      def weave(path)
        # This function works by moving through the selector path left-to-right,
        # building all possible prefixes simultaneously. These prefixes are
        # `befores`, while the remaining parenthesized suffixes is `afters`.
        befores = [[]]
        afters = path.dup

        until afters.empty?
          current = afters.shift.dup
          last_current = [current.pop]

          tempResult = []

          for before in befores do
            sub = subweave(before, current)
            if sub.nil?
              next
            end

            for seqs in sub do
              tempResult.push(seqs + last_current)
            end
          end

          befores = tempResult

        end

        return befores
      end
   */
  /*
      def weave(path)
        befores = [[]]
        afters = path.dup

        until afters.empty?
          current = afters.shift.dup

          last_current = [current.pop]


          tempResult = []

          for before in befores do
            sub = subweave(before, current)

            if sub.nil?
              next []
            end


            for seqs in sub do
              toPush = seqs + last_current

              tempResult.push(seqs + last_current)
            end

          end

          befores = tempResult

        end

        return befores
      end
  */
  Node Extend::weave(Node& path) {

    DEBUG_PRINTLN(WEAVE, "WEAVE: " << path)

    Node befores = Node::createCollection();
    befores.collection()->push_back(Node::createCollection());

    Node afters = Node::createCollection();
    afters.plus(path);

    while (!afters.collection()->empty()) {
      Node current = afters.collection()->front().klone();
      afters.collection()->pop_front();
      DEBUG_PRINTLN(WEAVE, "CURRENT: " << current)
      if (current.collection()->size() == 0) continue;

      Node last_current = Node::createCollection();
      last_current.collection()->push_back(current.collection()->back());
      current.collection()->pop_back();
      DEBUG_PRINTLN(WEAVE, "CURRENT POST POP: " << current)
      DEBUG_PRINTLN(WEAVE, "LAST CURRENT: " << last_current)

      Node tempResult = Node::createCollection();

      for (NodeDeque::iterator beforesIter = befores.collection()->begin(), beforesEndIter = befores.collection()->end(); beforesIter != beforesEndIter; beforesIter++) {
        Node& before = *beforesIter;

        Node sub = subweave(before, current);

        DEBUG_PRINTLN(WEAVE, "SUB: " << sub)

        if (sub.isNil()) {
          return Node::createCollection();
        }

        for (NodeDeque::iterator subIter = sub.collection()->begin(), subEndIter = sub.collection()->end(); subIter != subEndIter; subIter++) {
          Node& seqs = *subIter;

          Node toPush = Node::createCollection();
          toPush.plus(seqs);
          toPush.plus(last_current);

          // move line feed from inner to outer selector (very hacky indeed)
          if (last_current.collection() && last_current.collection()->front().selector()) {
            toPush.got_line_feed = last_current.collection()->front().got_line_feed;
            last_current.collection()->front().selector()->has_line_feed(false);
            last_current.collection()->front().got_line_feed = false;
          }

          tempResult.collection()->push_back(toPush);

        }
      }

      befores = tempResult;

    }

    return befores;
  }



  /*
   This is the equivalent of ruby's SimpleSequence.do_extend.

    // TODO: I think I have some modified ruby code to put here. Check.
  */
  /*
   ISSUES:
   - Previous TODO: Do we need to group the results by extender?
   - What does subject do in?: next unless unified = seq.members.last.unify(self_without_sel, subject?)
   - IMPROVEMENT: The search for uniqueness at the end is not ideal since it's has to loop over everything...
   - IMPROVEMENT: Check if the final search for uniqueness is doing anything that extendComplexSelector isn't already doing...
   */
  template<typename KeyType>
  class GroupByToAFunctor {
  public:
    KeyType operator()(SubSetMapPair& extPair) const {
      Complex_Selector_Obj pSelector = extPair.first;
      return pSelector;
    }
  };
  Node Extend::extendCompoundSelector(Compound_Selector_Ptr pSelector, CompoundSelectorSet& seen, bool isReplace) {

    /* this turned out to be too much overhead
       probably due to holding a "Node" object
    // check if we already extended this selector
    // we can do this since subset_map is "static"
    auto memoized = memoizeCompound.find(pSelector);
    if (memoized != memoizeCompound.end()) {
      return memoized->second.klone();
    }
    */

    DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pSelector, "EXTEND COMPOUND: "))
    // TODO: Ruby has another loop here to skip certain members?

    // let RESULTS be an empty list of complex selectors
    Node results = Node::createCollection();
    // extendedSelectors.got_line_feed = true;

    SubSetMapPairs entries = subset_map.get_v(pSelector);

    GroupByToAFunctor<Complex_Selector_Obj> extPairKeyFunctor;
    SubSetMapResults arr;
    group_by_to_a(entries, extPairKeyFunctor, arr);

    SubSetMapLookups holder;

    // for each (EXTENDER, TARGET) in MAP.get(COMPOUND):
    for (SubSetMapResult& groupedPair : arr) {

      Complex_Selector_Obj seq = groupedPair.first;
      SubSetMapPairs& group = groupedPair.second;

      DEBUG_EXEC(EXTEND_COMPOUND, printComplexSelector(seq, "SEQ: "))

      Compound_Selector_Obj pSels = SASS_MEMORY_NEW(Compound_Selector, pSelector->pstate());
      for (SubSetMapPair& pair : group) {
        pair.second->extended(true);
        pSels->concat(pair.second);
      }

      DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pSels, "SELS: "))

      // The selector up to where the @extend is (ie, the thing to merge)
      Complex_Selector_Ptr pExtComplexSelector = seq;

      // TODO: This can return a Compound_Selector with no elements. Should that just be returning NULL?
      // RUBY: self_without_sel = Sass::Util.array_minus(members, sels)
      Compound_Selector_Obj pSelectorWithoutExtendSelectors = pSelector->minus(pSels);

      DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pSelector, "MEMBERS: "))
      DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pSelectorWithoutExtendSelectors, "SELF_WO_SEL: "))

      Compound_Selector_Obj pInnermostCompoundSelector = pExtComplexSelector->last()->head();

      if (!pInnermostCompoundSelector) {
        pInnermostCompoundSelector = SASS_MEMORY_NEW(Compound_Selector, pSelector->pstate());
      }
      Compound_Selector_Obj pUnifiedSelector = pInnermostCompoundSelector->unify_with(pSelectorWithoutExtendSelectors);

      DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pInnermostCompoundSelector, "LHS: "))
      DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pSelectorWithoutExtendSelectors, "RHS: "))
      DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pUnifiedSelector, "UNIFIED: "))

      // RUBY: next unless unified
      if (!pUnifiedSelector || pUnifiedSelector->length() == 0) {
        continue;
      }

      // TODO: implement the parent directive match (if necessary based on test failures)
      // next if group.map {|e, _| check_directives_match!(e, parent_directives)}.none?

      // TODO: This seems a little fishy to me. See if it causes any problems. From the ruby, we should be able to just
      // get rid of the last Compound_Selector and replace it with this one. I think the reason this code is more
      // complex is that Complex_Selector contains a combinator, but in ruby combinators have already been filtered
      // out and aren't operated on.
      Complex_Selector_Obj pNewSelector = SASS_MEMORY_CLONE(pExtComplexSelector); // ->first();

      Complex_Selector_Obj pNewInnerMost = SASS_MEMORY_NEW(Complex_Selector, pSelector->pstate(), Complex_Selector::ANCESTOR_OF, pUnifiedSelector, {});

      Complex_Selector::Combinator combinator = pNewSelector->clear_innermost();
      pNewSelector->set_innermost(pNewInnerMost, combinator);

#ifdef DEBUG
      ComplexSelectorSet debugSet;
      debugSet = pNewSelector->sources();
      if (debugSet.size() > 0) {
        throw std::runtime_error("The new selector should start with no sources. Something needs to be cloned to fix this.");
      }
      debugSet = pExtComplexSelector->sources();
      if (debugSet.size() > 0) {
        throw std::runtime_error("The extension selector from our subset map should not have sources. These will bleed to the new selector. Something needs to be cloned to fix this.");
      }
#endif


      // if (pSelector && pSelector->has_line_feed()) pNewInnerMost->has_line_feed(true);
      // Set the sources on our new Complex_Selector to the sources of this simple sequence plus the thing we're extending.
      DEBUG_PRINTLN(EXTEND_COMPOUND, "SOURCES SETTING ON NEW SEQ: " << complexSelectorToNode(pNewSelector))

      DEBUG_EXEC(EXTEND_COMPOUND, ComplexSelectorSet oldSet = pNewSelector->sources(); printSourcesSet(oldSet, "SOURCES NEW SEQ BEGIN: "))

      // I actually want to create a copy here (performance!)
      ComplexSelectorSet newSourcesSet = pSelector->sources(); // XXX
      DEBUG_EXEC(EXTEND_COMPOUND, printSourcesSet(newSourcesSet, "SOURCES THIS EXTEND: "))

      newSourcesSet.insert(pExtComplexSelector);
      DEBUG_EXEC(EXTEND_COMPOUND, printSourcesSet(newSourcesSet, "SOURCES WITH NEW SOURCE: "))

      // RUBY: new_seq.add_sources!(sources + [seq])
      pNewSelector->addSources(newSourcesSet);

      DEBUG_EXEC(EXTEND_COMPOUND, ComplexSelectorSet newSet = pNewSelector->sources(); printSourcesSet(newSet, "SOURCES ON NEW SELECTOR AFTER ADD: "))
      DEBUG_EXEC(EXTEND_COMPOUND, printSourcesSet(pSelector->sources(), "SOURCES THIS EXTEND WHICH SHOULD BE SAME STILL: "))


      if (pSels->has_line_feed()) pNewSelector->has_line_feed(true);

      holder.push_back(std::make_pair(pSels, pNewSelector));
    }


    for (SubSetMapLookup& pair : holder) {

      Compound_Selector_Obj pSels = pair.first;
      Complex_Selector_Obj pNewSelector = pair.second;


      // RUBY??: next [] if seen.include?(sels)
      if (seen.find(pSels) != seen.end()) {
        continue;
      }


      CompoundSelectorSet recurseSeen(seen);
      recurseSeen.insert(pSels);


      DEBUG_PRINTLN(EXTEND_COMPOUND, "RECURSING DO EXTEND: " << complexSelectorToNode(pNewSelector))
      Node recurseExtendedSelectors = extendComplexSelector(pNewSelector, recurseSeen, isReplace, false); // !:isOriginal

      DEBUG_PRINTLN(EXTEND_COMPOUND, "RECURSING DO EXTEND RETURN: " << recurseExtendedSelectors)

      for (NodeDeque::iterator iterator = recurseExtendedSelectors.collection()->begin(), endIterator = recurseExtendedSelectors.collection()->end();
           iterator != endIterator; ++iterator) {
        Node newSelector = *iterator;

//        DEBUG_PRINTLN(EXTEND_COMPOUND, "EXTENDED AT THIS POINT: " << results)
//        DEBUG_PRINTLN(EXTEND_COMPOUND, "SELECTOR EXISTS ALREADY: " << newSelector << " " << results.contains(newSelector, false /*simpleSelectorOrderDependent*/));

        if (!results.contains(newSelector)) {
//          DEBUG_PRINTLN(EXTEND_COMPOUND, "ADDING NEW SELECTOR")
          results.collection()->push_back(newSelector);
        }
      }
    }

    DEBUG_EXEC(EXTEND_COMPOUND, printCompoundSelector(pSelector, "EXTEND COMPOUND END: "))

    // this turned out to be too much overhead
    // memory results in a map table - since extending is very expensive
    // memoizeCompound.insert(std::pair<Compound_Selector_Obj, Node>(pSelector, results));

    return results;
  }


  // check if selector has something to be extended by subset_map
  bool Extend::complexSelectorHasExtension(Complex_Selector_Ptr selector, CompoundSelectorSet& seen) {

    bool hasExtension = false;

    Complex_Selector_Obj pIter = selector;

    while (!hasExtension && pIter) {
      Compound_Selector_Obj pHead = pIter->head();

      if (pHead) {
        SubSetMapPairs entries = subset_map.get_v(pHead);
        for (SubSetMapPair ext : entries) {
          // check if both selectors have the same media block parent
          // if (ext.first->media_block() == pComplexSelector->media_block()) continue;
          if (ext.second->media_block() == 0) continue;
          if (pHead->media_block() &&
              ext.second->media_block()->media_queries() &&
              pHead->media_block()->media_queries()
          ) {
            std::string query_left(ext.second->media_block()->media_queries()->to_string());
            std::string query_right(pHead->media_block()->media_queries()->to_string());
            if (query_left == query_right) continue;
          }

          // fail if one goes across media block boundaries
          std::stringstream err;
          std::string cwd(Sass::File::get_cwd());
          ParserState pstate(ext.second->pstate());
          std::string rel_path(Sass::File::abs2rel(pstate.path, cwd, cwd));
          err << "You may not @extend an outer selector from within @media.\n";
          err << "You may only @extend selectors within the same directive.\n";
          err << "From \"@extend " << ext.second->to_string() << "\"";
          err << " on line " << pstate.line+1 << " of " << rel_path << "\n";
          error(err.str(), selector->pstate(), eval->exp.traces);
        }
        if (entries.size() > 0) hasExtension = true;
      }

      pIter = pIter->tail();
    }

    return hasExtension;
  }


  /*
   This is the equivalent of ruby's Sequence.do_extend.

   // TODO: I think I have some modified ruby code to put here. Check.
   */
  /*
   ISSUES:
   - check to automatically include combinators doesn't transfer over to libsass' data model where
     the combinator and compound selector are one unit
     next [[sseq_or_op]] unless sseq_or_op.is_a?(SimpleSequence)
   */
  Node Extend::extendComplexSelector(Complex_Selector_Ptr selector, CompoundSelectorSet& seen, bool isReplace, bool isOriginal) {

    // check if we already extended this selector
    // we can do this since subset_map is "static"
    auto memoized = memoizeComplex.find(selector);
    if (memoized != memoizeComplex.end()) {
      return memoized->second;
    }

    // convert the input selector to extend node format
    Node complexSelector = complexSelectorToNode(selector);
    DEBUG_PRINTLN(EXTEND_COMPLEX, "EXTEND COMPLEX: " << complexSelector)

    // let CHOICES be an empty list of selector-lists
    // create new collection to hold the results
    Node choices = Node::createCollection();

    // for each compound selector COMPOUND in COMPLEX:
    for (Node& sseqOrOp : *complexSelector.collection()) {

      DEBUG_PRINTLN(EXTEND_COMPLEX, "LOOP: " << sseqOrOp)

      // If it's not a selector (meaning it's a combinator), just include it automatically
      // RUBY: next [[sseq_or_op]] unless sseq_or_op.is_a?(SimpleSequence)
      if (!sseqOrOp.isSelector()) {
        // Wrap our Combinator in two collections to match ruby. This is essentially making a collection Node
        // with one collection child. The collection child represents a Complex_Selector that is only a combinator.
        Node outer = Node::createCollection();
        Node inner = Node::createCollection();
        outer.collection()->push_back(inner);
        inner.collection()->push_back(sseqOrOp);
        choices.collection()->push_back(outer);
        continue;
      }

      // verified now that node is a valid selector
      Complex_Selector_Obj sseqSel = sseqOrOp.selector();
      Compound_Selector_Obj sseqHead = sseqSel->head();

      // let EXTENDED be extend_compound(COMPOUND, SEEN)
      // extend the compound selector against the given subset_map
      // RUBY: extended = sseq_or_op.do_extend(extends, parent_directives, replace, seen)
      Node extended = extendCompoundSelector(sseqHead, seen, isReplace); // slow(17%)!
      if (sseqOrOp.got_line_feed) extended.got_line_feed = true;
      DEBUG_PRINTLN(EXTEND_COMPLEX, "EXTENDED: " << extended)

      // Prepend the Compound_Selector based on the choices logic; choices seems to be extend but with a ruby
      // Array instead of a Sequence due to the member mapping: choices = extended.map {|seq| seq.members}
      // RUBY: extended.first.add_sources!([self]) if original && !has_placeholder?
      if (isOriginal && !selector->has_placeholder()) {
        ComplexSelectorSet srcset;
        srcset.insert(selector);
        sseqSel->addSources(srcset);
        // DEBUG_PRINTLN(EXTEND_COMPLEX, "ADD SOURCES: " << *pComplexSelector)
      }

      bool isSuperselector = false;
      // if no complex selector in EXTENDED is a superselector of COMPOUND:
      for (Node& childNode : *extended.collection()) {
        Complex_Selector_Obj pExtensionSelector = nodeToComplexSelector(childNode);
        if (pExtensionSelector->is_superselector_of(sseqSel)) {
          isSuperselector = true;
          break;
        }
      }

      if (!isSuperselector) {
        // add a complex selector composed only of COMPOUND to EXTENDED
        if (sseqOrOp.got_line_feed) sseqSel->has_line_feed(sseqOrOp.got_line_feed);
        extended.collection()->push_front(complexSelectorToNode(sseqSel));
      }

      DEBUG_PRINTLN(EXTEND_COMPLEX, "CHOICES UNSHIFTED: " << extended)

      // add EXTENDED to CHOICES
      // Aggregate our current extensions
      choices.collection()->push_back(extended);
    }


    DEBUG_PRINTLN(EXTEND_COMPLEX, "EXTENDED NOT EXPANDED: " << choices)



    // Ruby Equivalent: paths
    Node paths = Sass::paths(choices);

    DEBUG_PRINTLN(EXTEND_COMPLEX, "PATHS: " << paths)

    // let WEAVES be an empty list of selector lists
    Node weaves = Node::createCollection();

    // for each list of complex selectors PATH in paths(CHOICES):
    for (Node& path : *paths.collection()) {
      // add weave(PATH) to WEAVES
      Node weaved = weave(path); // slow(12%)!
      weaved.got_line_feed = path.got_line_feed;
      weaves.collection()->push_back(weaved);
    }

    DEBUG_PRINTLN(EXTEND_COMPLEX, "WEAVES: " << weaves)

    // Ruby Equivalent: trim
    Node trimmed(trim(weaves, isReplace)); // slow(19%)!

    DEBUG_PRINTLN(EXTEND_COMPLEX, "TRIMMED: " << trimmed)

    // Ruby Equivalent: flatten
    Node flattened(flatten(trimmed, 1));

    DEBUG_PRINTLN(EXTEND_COMPLEX, "FLATTENED: " << flattened)

    // memory results in a map table - since extending is very expensive
    memoizeComplex.insert(std::pair<Complex_Selector_Obj, Node>(selector, flattened));

    // return trim(WEAVES)
    return flattened;
  }



  /*
   This is the equivalent of ruby's CommaSequence.do_extend.
  */
  // We get a selector list with has something to extend and a subset_map with
  // all extenders. Pick the ones that match our selectors in the list.
  Selector_List_Ptr Extend::extendSelectorList(Selector_List_Obj pSelectorList, bool isReplace, bool& extendedSomething, CompoundSelectorSet& seen) {

    Selector_List_Obj pNewSelectors = SASS_MEMORY_NEW(Selector_List, pSelectorList->pstate(), pSelectorList->length());

    // check if we already extended this selector
    // we can do this since subset_map is "static"
    auto memoized = memoizeList.find(pSelectorList);
    if (memoized != memoizeList.end()) {
      extendedSomething = true;
      return memoized->second;
    }

    extendedSomething = false;
    // process each comlplex selector in the selector list.
    // Find the ones that can be extended by given subset_map.
    for (size_t index = 0, length = pSelectorList->length(); index < length; index++) {
      Complex_Selector_Obj pSelector = (*pSelectorList)[index];

      // ruby sass seems to keep a list of things that have extensions and then only extend those. We don't currently do that.
      // Since it's not that expensive to check if an extension exists in the subset map and since it can be relatively expensive to
      // run through the extend code (which does a data model transformation), check if there is anything to extend before doing
      // the extend. We might be able to optimize extendComplexSelector, but this approach keeps us closer to ruby sass (which helps
      // when debugging).
      if (!complexSelectorHasExtension(pSelector, seen)) {
        pNewSelectors->append(pSelector);
        continue;
      }

      // complexSelectorHasExtension was true!
      extendedSomething = true;

      // now do the actual extension of the complex selector
      Node extendedSelectors = extendComplexSelector(pSelector, seen, isReplace, true);

      if (!pSelector->has_placeholder()) {
        Node nSelector(complexSelectorToNode(pSelector));
        if (!extendedSelectors.contains(nSelector)) {
          pNewSelectors->append(pSelector);
          continue;
        }
      }

      bool doReplace = isReplace;
      for (Node& childNode : *extendedSelectors.collection()) {
        // When it is a replace, skip the first one, unless there is only one
        if(doReplace && extendedSelectors.collection()->size() > 1 ) {
          doReplace = false;
          continue;
        }
        pNewSelectors->append(nodeToComplexSelector(childNode));
      }
    }

    Remove_Placeholders remove_placeholders;
    // it seems that we have to remove the place holders early here
    // normally we do this as the very last step (compare to ruby sass)
    pNewSelectors = remove_placeholders.remove_placeholders(pNewSelectors);

    // unwrap all wrapped selectors with inner lists
    for (Complex_Selector_Obj cur : pNewSelectors->elements()) {
      // process tails
      while (cur) {
        // process header
        if (cur->head() && seen.find(cur->head()) == seen.end()) {
          CompoundSelectorSet recseen(seen);
          recseen.insert(cur->head());
          // create a copy since we add multiple items if stuff get unwrapped
          Compound_Selector_Obj cpy_head = SASS_MEMORY_NEW(Compound_Selector, cur->pstate());
          for (Simple_Selector_Obj hs : *cur->head()) {
            if (Wrapped_Selector_Obj ws = Cast<Wrapped_Selector>(hs)) {
              ws->selector(SASS_MEMORY_CLONE(ws->selector()));
              if (Selector_List_Obj sl = Cast<Selector_List>(ws->selector())) {
                // special case for ruby ass
                if (sl->empty()) {
                  // this seems inconsistent but it is how ruby sass seems to remove parentheses
                  cpy_head->append(SASS_MEMORY_NEW(Type_Selector, hs->pstate(), ws->name()));
                }
                // has wrapped not selectors
                else if (ws->name() == ":not") {
                  // extend the inner list of wrapped selector
                  bool extended = false;
                  Selector_List_Obj ext_sl = extendSelectorList(sl, false, extended, recseen);
                  for (size_t i = 0; i < ext_sl->length(); i += 1) {
                    if (Complex_Selector_Obj ext_cs = ext_sl->at(i)) {
                      // create clones for wrapped selector and the inner list
                      Wrapped_Selector_Obj cpy_ws = SASS_MEMORY_COPY(ws);
                      Selector_List_Obj cpy_ws_sl = SASS_MEMORY_NEW(Selector_List, sl->pstate());
                      // remove parent selectors from inner selector
                      Compound_Selector_Obj ext_head;
                      if (ext_cs->first()) ext_head = ext_cs->first()->head();
                      if (ext_head && ext_head && ext_head->length() > 0) {
                        cpy_ws_sl->append(ext_cs->mutable_first());
                      }
                      // assign list to clone
                      cpy_ws->selector(cpy_ws_sl);
                      // append the clone
                      cpy_head->append(cpy_ws);
                    }
                  }
                  if (eval && extended) {
                    eval->exp.selector_stack.push_back(pNewSelectors);
                    cpy_head->perform(eval);
                    eval->exp.selector_stack.pop_back();
                  }
                }
                // has wrapped selectors
                else {
                  Wrapped_Selector_Obj cpy_ws = SASS_MEMORY_COPY(ws);
                  Selector_List_Obj ext_sl = extendSelectorList(sl, recseen);
                  cpy_ws->selector(ext_sl);
                  cpy_head->append(cpy_ws);
                }
              } else {
                cpy_head->append(hs);
              }
            } else {
              cpy_head->append(hs);
            }
          }
          // replace header
          cur->head(cpy_head);
        }
        // process tail
        cur = cur->tail();
      }
    }

    // memory results in a map table - since extending is very expensive
    memoizeList.insert(std::pair<Selector_List_Obj, Selector_List_Obj>(pSelectorList, pNewSelectors));

    return pNewSelectors.detach();

  }


  bool shouldExtendBlock(Block_Obj b) {

    // If a block is empty, there's no reason to extend it since any rules placed on this block
    // won't have any output. The main benefit of this is for structures like:
    //
    //    .a {
    //      .b {
    //        x: y;
    //      }
    //    }
    //
    // We end up visiting two rulesets (one with the selector .a and the other with the selector .a .b).
    // In this case, we don't want to try to pull rules onto .a since they won't get output anyway since
    // there are no child statements. However .a .b should have extensions applied.

    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj stm = b->at(i);

      if (Cast<Ruleset>(stm)) {
        // Do nothing. This doesn't count as a statement that causes extension since we'll
        // iterate over this rule set in a future visit and try to extend it.
      }
      else {
        return true;
      }
    }

    return false;

  }


  // Extend a ruleset by extending the selectors and updating them on the ruleset. The block's rules don't need to change.
  // Every Ruleset in the whole tree is calling this function. We decide if there
  // was is @extend that matches our selector. If we find one, we will go further
  // and call the extend magic for our selector. The subset_map contains all blocks
  // where @extend was found. Pick the ones that match our selector!
  void Extend::extendObjectWithSelectorAndBlock(Ruleset_Ptr pObject) {

    DEBUG_PRINTLN(EXTEND_OBJECT, "FOUND SELECTOR: " << Cast<Selector_List>(pObject->selector())->to_string())

    // Ruby sass seems to filter nodes that don't have any content well before we get here.
    // I'm not sure the repercussions of doing so, so for now, let's just not extend things
    // that won't be output later. Profiling shows this may us 0.2% or so.
    if (!shouldExtendBlock(pObject->block())) {
      DEBUG_PRINTLN(EXTEND_OBJECT, "RETURNING WITHOUT EXTEND ATTEMPT")
      return;
    }

    bool extendedSomething = false;

    CompoundSelectorSet seen;
    Selector_List_Obj pNewSelectorList = extendSelectorList(pObject->selector(), false, extendedSomething, seen);

    if (extendedSomething && pNewSelectorList) {
      DEBUG_PRINTLN(EXTEND_OBJECT, "EXTEND ORIGINAL SELECTORS: " << pObject->selector()->to_string())
      DEBUG_PRINTLN(EXTEND_OBJECT, "EXTEND SETTING NEW SELECTORS: " << pNewSelectorList->to_string())
      pNewSelectorList->remove_parent_selectors();
      pObject->selector(pNewSelectorList);
    } else {
      DEBUG_PRINTLN(EXTEND_OBJECT, "EXTEND DID NOT TRY TO EXTEND ANYTHING")
    }
  }

  Extend::Extend(Subset_Map& ssm)
  : subset_map(ssm), eval(NULL)
  { }

  void Extend::setEval(Eval& e) {
    eval = &e;
  }

  void Extend::operator()(Block_Ptr b)
  {
    for (size_t i = 0, L = b->length(); i < L; ++i) {
      Statement_Obj stm = b->at(i);
      stm->perform(this);
    }
    // do final check if everything was extended
    // we set `extended` flag on extended selectors
    if (b->is_root()) {
      // debug_subset_map(subset_map);
      for(auto const &it : subset_map.values()) {
        Complex_Selector_Ptr_Const sel = nullptr;
        Compound_Selector_Ptr_Const ext = nullptr;
        if (it.first) sel = it.first->first();
        if (it.second) ext = it.second;
        if (ext && (ext->extended() || ext->is_optional())) continue;
        std::string str_sel(sel ? sel->to_string({ NESTED, 5 }) : "NULL");
        std::string str_ext(ext ? ext->to_string({ NESTED, 5 }) : "NULL");
        // debug_ast(sel, "sel: ");
        // debug_ast(ext, "ext: ");
        error("\"" + str_sel + "\" failed to @extend \"" + str_ext + "\".\n"
              "The selector \"" + str_ext + "\" was not found.\n"
              "Use \"@extend " + str_ext + " !optional\" if the"
              " extend should be able to fail.", (ext ? ext->pstate() : NULL), eval->exp.traces);
      }
    }

  }

  void Extend::operator()(Ruleset_Ptr pRuleset)
  {
    extendObjectWithSelectorAndBlock( pRuleset );
    pRuleset->block()->perform(this);
  }

  void Extend::operator()(Supports_Block_Ptr pFeatureBlock)
  {
    pFeatureBlock->block()->perform(this);
  }

  void Extend::operator()(Media_Block_Ptr pMediaBlock)
  {
    pMediaBlock->block()->perform(this);
  }

  void Extend::operator()(Directive_Ptr a)
  {
    // Selector_List_Ptr ls = Cast<Selector_List>(a->selector());
    // selector_stack.push_back(ls);
    if (a->block()) a->block()->perform(this);
    // exp.selector_stack.pop_back();
  }
}

#if defined(_MSC_VER) && _MSC_VER >= 1800 && _MSC_VER < 1900 && defined(_M_X64)
#include <mutex>
#endif

namespace Sass {

  double round(double val, size_t precision)
  {
    // Disable FMA3-optimized implementation when compiling with VS2013 for x64 targets
    // See https://github.com/sass/node-sass/issues/1854 for details
    // FIXME: Remove this workaround when we switch to VS2015+
    #if defined(_MSC_VER) && _MSC_VER >= 1800 && _MSC_VER < 1900 && defined(_M_X64)
      static std::once_flag flag;
      std::call_once(flag, []() { _set_FMA3_enable(0); });
    #endif

    // https://github.com/sass/sass/commit/4e3e1d5684cc29073a507578fc977434ff488c93
    if (fmod(val, 1) - 0.5 > - std::pow(0.1, precision + 1)) return std::ceil(val);
    else if (fmod(val, 1) - 0.5 > std::pow(0.1, precision)) return std::floor(val);
    // work around some compiler issue
    // cygwin has it not defined in std
    using namespace std;
    return ::round(val);
  }

  /* Locale unspecific atof function. */
  double sass_strtod(const char *str)
  {
    char separator = *(localeconv()->decimal_point);
    if(separator != '.'){
      // The current locale specifies another
      // separator. convert the separator to the
      // one understood by the locale if needed
      const char *found = strchr(str, '.');
      if(found != NULL){
        // substitution is required. perform the substitution on a copy
        // of the string. This is slower but it is thread safe.
        char *copy = sass_copy_c_string(str);
        *(copy + (found - str)) = separator;
        double res = strtod(copy, NULL);
        free(copy);
        return res;
      }
    }

    return strtod(str, NULL);
  }

  // helper for safe access to c_ctx
  const char* safe_str (const char* str, const char* alt) {
    return str == NULL ? alt : str;
  }

  void free_string_array(char ** arr) {
    if(!arr)
        return;

    char **it = arr;
    while (it && (*it)) {
      free(*it);
      ++it;
    }

    free(arr);
  }

  char **copy_strings(const std::vector<std::string>& strings, char*** array, int skip) {
    int num = static_cast<int>(strings.size()) - skip;
    char** arr = (char**) calloc(num + 1, sizeof(char*));
    if (arr == 0)
      return *array = (char **)NULL;

    for(int i = 0; i < num; i++) {
      arr[i] = (char*) malloc(sizeof(char) * (strings[i + skip].size() + 1));
      if (arr[i] == 0) {
        free_string_array(arr);
        return *array = (char **)NULL;
      }
      std::copy(strings[i + skip].begin(), strings[i + skip].end(), arr[i]);
      arr[i][strings[i + skip].size()] = '\0';
    }

    arr[num] = 0;
    return *array = arr;
  }

  // read css string (handle multiline DELIM)
  std::string read_css_string(const std::string& str, bool css)
  {
    if (!css) return str;
    std::string out("");
    bool esc = false;
    for (auto i : str) {
      if (i == '\\') {
        esc = ! esc;
      } else if (esc && i == '\r') {
        continue;
      } else if (esc && i == '\n') {
        out.resize (out.size () - 1);
        esc = false;
        continue;
      } else {
        esc = false;
      }
      out.push_back(i);
    }
    // happens when parsing does not correctly skip
    // over escaped sequences for ie. interpolations
    // one example: foo\#{interpolate}
    // if (esc) out += '\\';
    return out;
  }

  // double escape all escape sequences
  // keep unescaped quotes and backslashes
  std::string evacuate_escapes(const std::string& str)
  {
    std::string out("");
    bool esc = false;
    for (auto i : str) {
      if (i == '\\' && !esc) {
        out += '\\';
        out += '\\';
        esc = true;
      } else if (esc && i == '"') {
        out += '\\';
        out += i;
        esc = false;
      } else if (esc && i == '\'') {
        out += '\\';
        out += i;
        esc = false;
      } else if (esc && i == '\\') {
        out += '\\';
        out += i;
        esc = false;
      } else {
        esc = false;
        out += i;
      }
    }
    // happens when parsing does not correctly skip
    // over escaped sequences for ie. interpolations
    // one example: foo\#{interpolate}
    // if (esc) out += '\\';
    return out;
  }

  // bell characters are replaced with spaces
  void newline_to_space(std::string& str)
  {
    std::replace(str.begin(), str.end(), '\n', ' ');
  }

  // bell characters are replaced with spaces
  // also eats spaces after line-feeds (ltrim)
  std::string string_to_output(const std::string& str)
  {
    std::string out("");
    bool lf = false;
    for (auto i : str) {
      if (i == '\n') {
        out += ' ';
        lf = true;
      } else if (!(lf && isspace(i))) {
        out += i;
        lf = false;
      }
    }
    return out;
  }

  std::string escape_string(const std::string& str)
  {
    std::string out("");
    for (auto i : str) {
      if (i == '\n') {
        out += "\\n";
      } else if (i == '\r') {
        out += "\\r";
      } else if (i == '\t') {
        out += "\\t";
      } else {
        out += i;
      }
    }
    return out;
  }

  std::string comment_to_string(const std::string& text)
  {
    std::string str = "";
    size_t has = 0;
    char prev = 0;
    bool clean = false;
    for (auto i : text) {
      if (clean) {
        if (i == '\n') { has = 0; }
        else if (i == '\r') { has = 0; }
        else if (i == '\t') { ++ has; }
        else if (i == ' ') { ++ has; }
        else if (i == '*') {}
        else {
          clean = false;
          str += ' ';
          if (prev == '*' && i == '/') str += "*/";
          else str += i;
        }
      } else if (i == '\n') {
        clean = true;
      } else if (i == '\r') {
        clean = true;
      } else {
        str += i;
      }
      prev = i;
    }
    if (has) return str;
    else return text;
  }

  // find best quote_mark by detecting if the string contains any single
  // or double quotes. When a single quote is found, we not we want a double
  // quote as quote_mark. Otherwise we check if the string cotains any double
  // quotes, which will trigger the use of single quotes as best quote_mark.
  char detect_best_quotemark(const char* s, char qm)
  {
    // ensure valid fallback quote_mark
    char quote_mark = qm && qm != '*' ? qm : '"';
    while (*s) {
      // force double quotes as soon
      // as one single quote is found
      if (*s == '\'') { return '"'; }
      // a single does not force quote_mark
      // maybe we see a double quote later
      else if (*s == '"') { quote_mark = '\''; }
      ++ s;
    }
    return quote_mark;
  }

  std::string read_hex_escapes(const std::string& s)
  {

    std::string result;
    bool skipped = false;

    for (size_t i = 0, L = s.length(); i < L; ++i) {

      // implement the same strange ruby sass behavior
      // an escape sequence can also mean a unicode char
      if (s[i] == '\\' && !skipped) {

        // remember
        skipped = true;

        // escape length
        size_t len = 1;

        // parse as many sequence chars as possible
        // ToDo: Check if ruby aborts after possible max
        while (i + len < L && s[i + len] && isxdigit(s[i + len])) ++ len;

        if (len > 1) {

          // convert the extracted hex string to code point value
          // ToDo: Maybe we could do this without creating a substring
          uint32_t cp = strtol(s.substr (i + 1, len - 1).c_str(), NULL, 16);

          if (s[i + len] == ' ') ++ len;

          // assert invalid code points
          if (cp == 0) cp = 0xFFFD;
          // replace bell character
          // if (cp == '\n') cp = 32;

          // use a very simple approach to convert via utf8 lib
          // maybe there is a more elegant way; maybe we shoud
          // convert the whole output from string to a stream!?
          // allocate memory for utf8 char and convert to utf8
          unsigned char u[5] = {0,0,0,0,0}; utf8::append(cp, u);
          for(size_t m = 0; m < 5 && u[m]; m++) result.push_back(u[m]);

          // skip some more chars?
          i += len - 1; skipped = false;

        }

        else {

          skipped = false;

          result.push_back(s[i]);

        }

      }

      else {

        result.push_back(s[i]);

      }

    }

    return result;

  }

  std::string unquote(const std::string& s, char* qd, bool keep_utf8_sequences, bool strict)
  {

    // not enough room for quotes
    // no possibility to unquote
    if (s.length() < 2) return s;

    char q;
    bool skipped = false;

    // this is no guarantee that the unquoting will work
    // what about whitespace before/after the quote_mark?
    if      (*s.begin() == '"'  && *s.rbegin() == '"')  q = '"';
    else if (*s.begin() == '\'' && *s.rbegin() == '\'') q = '\'';
    else                                                return s;

    std::string unq;
    unq.reserve(s.length()-2);

    for (size_t i = 1, L = s.length() - 1; i < L; ++i) {

      // implement the same strange ruby sass behavior
      // an escape sequence can also mean a unicode char
      if (s[i] == '\\' && !skipped) {
        // remember
        skipped = true;

        // skip it
        // ++ i;

        // if (i == L) break;

        // escape length
        size_t len = 1;

        // parse as many sequence chars as possible
        // ToDo: Check if ruby aborts after possible max
        while (i + len < L && s[i + len] && isxdigit(s[i + len])) ++ len;

        // hex string?
        if (keep_utf8_sequences) {
          unq.push_back(s[i]);
        } else if (len > 1) {

          // convert the extracted hex string to code point value
          // ToDo: Maybe we could do this without creating a substring
          uint32_t cp = strtol(s.substr (i + 1, len - 1).c_str(), NULL, 16);

          if (s[i + len] == ' ') ++ len;

          // assert invalid code points
          if (cp == 0) cp = 0xFFFD;
          // replace bell character
          // if (cp == '\n') cp = 32;

          // use a very simple approach to convert via utf8 lib
          // maybe there is a more elegant way; maybe we shoud
          // convert the whole output from string to a stream!?
          // allocate memory for utf8 char and convert to utf8
          unsigned char u[5] = {0,0,0,0,0}; utf8::append(cp, u);
          for(size_t m = 0; m < 5 && u[m]; m++) unq.push_back(u[m]);

          // skip some more chars?
          i += len - 1; skipped = false;

        }


      }
      // check for unexpected delimiter
      // be strict and throw error back
      // else if (!skipped && q == s[i]) {
      //   // don't be that strict
      //   return s;
      //   // this basically always means an internal error and not users fault
      //   error("Unescaped delimiter in string to unquote found. [" + s + "]", ParserState("[UNQUOTE]"));
      // }
      else {
        if (strict && !skipped) {
          if (s[i] == q) return s;
        }
        skipped = false;
        unq.push_back(s[i]);
      }

    }
    if (skipped) { return s; }
    if (qd) *qd = q;
    return unq;

  }

  std::string quote(const std::string& s, char q)
  {

    // autodetect with fallback to given quote
    q = detect_best_quotemark(s.c_str(), q);

    // return an empty quoted string
    if (s.empty()) return std::string(2, q ? q : '"');

    std::string quoted;
    quoted.reserve(s.length()+2);
    quoted.push_back(q);

    const char* it = s.c_str();
    const char* end = it + strlen(it) + 1;
    while (*it && it < end) {
      const char* now = it;

      if (*it == q) {
        quoted.push_back('\\');
      } else if (*it == '\\') {
        quoted.push_back('\\');
      }

      int cp = utf8::next(it, end);

      // in case of \r, check if the next in sequence
      // is \n and then advance the iterator and skip \r
      if (cp == '\r' && it < end && utf8::peek_next(it, end) == '\n') {
        cp = utf8::next(it, end);
      }

      if (cp == '\n') {
        quoted.push_back('\\');
        quoted.push_back('a');
        // we hope we can remove this flag once we figure out
        // why ruby sass has these different output behaviors
        // gsub(/\n(?![a-fA-F0-9\s])/, "\\a").gsub("\n", "\\a ")
        using namespace Prelexer;
        if (alternatives <
          Prelexer::char_range<'a', 'f'>,
          Prelexer::char_range<'A', 'F'>,
          Prelexer::char_range<'0', '9'>,
          space
        >(it) != NULL) {
          quoted.push_back(' ');
        }
      } else if (cp < 127) {
        quoted.push_back((char) cp);
      } else {
        while (now < it) {
          quoted.push_back(*now);
          ++ now;
        }
      }
    }

    quoted.push_back(q);
    return quoted;
  }

  bool is_hex_doublet(double n)
  {
    return n == 0x00 || n == 0x11 || n == 0x22 || n == 0x33 ||
           n == 0x44 || n == 0x55 || n == 0x66 || n == 0x77 ||
           n == 0x88 || n == 0x99 || n == 0xAA || n == 0xBB ||
           n == 0xCC || n == 0xDD || n == 0xEE || n == 0xFF ;
  }

  bool is_color_doublet(double r, double g, double b)
  {
    return is_hex_doublet(r) && is_hex_doublet(g) && is_hex_doublet(b);
  }

  bool peek_linefeed(const char* start)
  {
    using namespace Prelexer;
    using namespace Constants;
    return sequence <
             zero_plus <
               alternatives <
                 exactly <' '>,
                 exactly <'\t'>,
                 line_comment,
                 block_comment,
                 delimited_by <
                   slash_star,
                   star_slash,
                   false
                 >
               >
             >,
             re_linebreak
           >(start) != 0;
  }

  namespace Util {
    using std::string;

    std::string rtrim(const std::string &str) {
      std::string trimmed = str;
      size_t pos_ws = trimmed.find_last_not_of(" \t\n\v\f\r");
      if (pos_ws != std::string::npos)
      { trimmed.erase(pos_ws + 1); }
      else { trimmed.clear(); }
      return trimmed;
    }

    std::string normalize_underscores(const std::string& str) {
      std::string normalized = str;
      for(size_t i = 0, L = normalized.length(); i < L; ++i) {
        if(normalized[i] == '_') {
          normalized[i] = '-';
        }
      }
      return normalized;
    }

    std::string normalize_decimals(const std::string& str) {
      std::string prefix = "0";
      std::string normalized = str;

      return normalized[0] == '.' ? normalized.insert(0, prefix) : normalized;
    }

    bool isPrintable(Ruleset_Ptr r, Sass_Output_Style style) {
      if (r == NULL) {
        return false;
      }

      Block_Obj b = r->block();

      Selector_List_Ptr sl = Cast<Selector_List>(r->selector());
      bool hasSelectors = sl ? sl->length() > 0 : false;

      if (!hasSelectors) {
        return false;
      }

      bool hasDeclarations = false;
      bool hasPrintableChildBlocks = false;
      for (size_t i = 0, L = b->length(); i < L; ++i) {
        Statement_Obj stm = b->at(i);
        if (Cast<Directive>(stm)) {
          return true;
        } else if (Declaration_Ptr d = Cast<Declaration>(stm)) {
          return isPrintable(d, style);
        } else if (Has_Block_Ptr p = Cast<Has_Block>(stm)) {
          Block_Obj pChildBlock = p->block();
          if (isPrintable(pChildBlock, style)) {
            hasPrintableChildBlocks = true;
          }
        } else if (Comment_Ptr c = Cast<Comment>(stm)) {
          // keep for uncompressed
          if (style != COMPRESSED) {
            hasDeclarations = true;
          }
          // output style compressed
          if (c->is_important()) {
            hasDeclarations = c->is_important();
          }
        } else {
          hasDeclarations = true;
        }

        if (hasDeclarations || hasPrintableChildBlocks) {
          return true;
        }
      }

      return false;
    }

    bool isPrintable(String_Constant_Ptr s, Sass_Output_Style style)
    {
      return ! s->value().empty();
    }

    bool isPrintable(String_Quoted_Ptr s, Sass_Output_Style style)
    {
      return true;
    }

    bool isPrintable(Declaration_Ptr d, Sass_Output_Style style)
    {
      Expression_Obj val = d->value();
      if (String_Quoted_Obj sq = Cast<String_Quoted>(val)) return isPrintable(sq.ptr(), style);
      if (String_Constant_Obj sc = Cast<String_Constant>(val)) return isPrintable(sc.ptr(), style);
      return true;
    }

    bool isPrintable(Supports_Block_Ptr f, Sass_Output_Style style) {
      if (f == NULL) {
        return false;
      }

      Block_Obj b = f->block();

      bool hasDeclarations = false;
      bool hasPrintableChildBlocks = false;
      for (size_t i = 0, L = b->length(); i < L; ++i) {
        Statement_Obj stm = b->at(i);
        if (Cast<Declaration>(stm) || Cast<Directive>(stm)) {
          hasDeclarations = true;
        }
        else if (Has_Block_Ptr b = Cast<Has_Block>(stm)) {
          Block_Obj pChildBlock = b->block();
          if (!b->is_invisible()) {
            if (isPrintable(pChildBlock, style)) {
              hasPrintableChildBlocks = true;
            }
          }
        }

        if (hasDeclarations || hasPrintableChildBlocks) {
          return true;
        }
      }

      return false;
    }

    bool isPrintable(Media_Block_Ptr m, Sass_Output_Style style)
    {
      if (m == 0) return false;
      Block_Obj b = m->block();
      if (b == 0) return false;
      for (size_t i = 0, L = b->length(); i < L; ++i) {
        Statement_Obj stm = b->at(i);
        if (Cast<Directive>(stm)) return true;
        else if (Cast<Declaration>(stm)) return true;
        else if (Comment_Ptr c = Cast<Comment>(stm)) {
          if (isPrintable(c, style)) {
            return true;
          }
        }
        else if (Ruleset_Ptr r = Cast<Ruleset>(stm)) {
          if (isPrintable(r, style)) {
            return true;
          }
        }
        else if (Supports_Block_Ptr f = Cast<Supports_Block>(stm)) {
          if (isPrintable(f, style)) {
            return true;
          }
        }
        else if (Media_Block_Ptr mb = Cast<Media_Block>(stm)) {
          if (isPrintable(mb, style)) {
            return true;
          }
        }
        else if (Has_Block_Ptr b = Cast<Has_Block>(stm)) {
          if (isPrintable(b->block(), style)) {
            return true;
          }
        }
      }
      return false;
    }

    bool isPrintable(Comment_Ptr c, Sass_Output_Style style)
    {
      // keep for uncompressed
      if (style != COMPRESSED) {
        return true;
      }
      // output style compressed
      if (c->is_important()) {
        return true;
      }
      // not printable
      return false;
    };

    bool isPrintable(Block_Obj b, Sass_Output_Style style) {
      if (!b) {
        return false;
      }

      for (size_t i = 0, L = b->length(); i < L; ++i) {
        Statement_Obj stm = b->at(i);
        if (Cast<Declaration>(stm) || Cast<Directive>(stm)) {
          return true;
        }
        else if (Comment_Ptr c = Cast<Comment>(stm)) {
          if (isPrintable(c, style)) {
            return true;
          }
        }
        else if (Ruleset_Ptr r = Cast<Ruleset>(stm)) {
          if (isPrintable(r, style)) {
            return true;
          }
        }
        else if (Supports_Block_Ptr f = Cast<Supports_Block>(stm)) {
          if (isPrintable(f, style)) {
            return true;
          }
        }
        else if (Media_Block_Ptr m = Cast<Media_Block>(stm)) {
          if (isPrintable(m, style)) {
            return true;
          }
        }
        else if (Has_Block_Ptr b = Cast<Has_Block>(stm)) {
          if (isPrintable(b->block(), style)) {
            return true;
          }
        }
      }

      return false;
    }

    bool isAscii(const char chr) {
      return unsigned(chr) < 128;
    }

  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


#ifdef _WIN32
#include <windows.h>
#else
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <dlfcn.h>
#endif

namespace Sass {

  Plugins::Plugins(void) { }
  Plugins::~Plugins(void)
  {
    for (auto function : functions) {
      sass_delete_function(function);
    }
    for (auto importer : importers) {
      sass_delete_importer(importer);
    }
    for (auto header : headers) {
      sass_delete_importer(header);
    }
  }

  // check if plugin is compatible with this version
  // plugins may be linked static against libsass
  // we try to be compatible between major versions
  inline bool compatibility(const char* their_version)
  {
// const char* their_version = "3.1.2";
    // first check if anyone has an unknown version
    const char* our_version = libsass_version();
    if (!strcmp(their_version, "[na]")) return false;
    if (!strcmp(our_version, "[na]")) return false;

    // find the position of the second dot
    size_t pos = std::string(our_version).find('.', 0);
    if (pos != std::string::npos) pos = std::string(our_version).find('.', pos + 1);

    // if we do not have two dots we fallback to compare complete string
    if (pos == std::string::npos) { return strcmp(their_version, our_version) ? 0 : 1; }
    // otherwise only compare up to the second dot (major versions)
    else { return strncmp(their_version, our_version, pos) ? 0 : 1; }

  }

  // load one specific plugin
  bool Plugins::load_plugin (const std::string& path)
  {

    typedef const char* (*__plugin_version__)(void);
    typedef Sass_Function_List (*__plugin_load_fns__)(void);
    typedef Sass_Importer_List (*__plugin_load_imps__)(void);

    if (LOAD_LIB(plugin, path))
    {
      // try to load initial function to query libsass version suppor
      if (LOAD_LIB_FN(__plugin_version__, plugin_version, "libsass_get_version"))
      {
        // get the libsass version of the plugin
        if (!compatibility(plugin_version())) return false;
        // try to get import address for "libsass_load_functions"
        if (LOAD_LIB_FN(__plugin_load_fns__, plugin_load_functions, "libsass_load_functions"))
        {
          Sass_Function_List fns = plugin_load_functions(), _p = fns;
          while (fns && *fns) { functions.push_back(*fns); ++ fns; }
          sass_free_memory(_p); // only delete the container, items not yet
        }
        // try to get import address for "libsass_load_importers"
        if (LOAD_LIB_FN(__plugin_load_imps__, plugin_load_importers, "libsass_load_importers"))
        {
          Sass_Importer_List imps = plugin_load_importers(), _p = imps;
          while (imps && *imps) { importers.push_back(*imps); ++ imps; }
          sass_free_memory(_p); // only delete the container, items not yet
        }
        // try to get import address for "libsass_load_headers"
        if (LOAD_LIB_FN(__plugin_load_imps__, plugin_load_headers, "libsass_load_headers"))
        {
          Sass_Importer_List imps = plugin_load_headers(), _p = imps;
          while (imps && *imps) { headers.push_back(*imps); ++ imps; }
          sass_free_memory(_p); // only delete the container, items not yet
        }
        // success
        return true;
      }
      else
      {
        // print debug message to stderr (should not happen)
        std::cerr << "failed loading 'libsass_support' in <" << path << ">" << std::endl;
        if (const char* dlsym_error = dlerror()) std::cerr << dlsym_error << std::endl;
        CLOSE_LIB(plugin);
      }
    }
    else
    {
      // print debug message to stderr (should not happen)
      std::cerr << "failed loading plugin <" << path << ">" << std::endl;
      if (const char* dlopen_error = dlerror()) std::cerr << dlopen_error << std::endl;
    }

    return false;

  }

  size_t Plugins::load_plugins(const std::string& path)
  {

    // count plugins
    size_t loaded = 0;

    #ifdef _WIN32

      try
      {

        // use wchar (utf16)
        WIN32_FIND_DATAW data;
        // trailing slash is guaranteed
        std::string globsrch(path + "*.dll");
        // convert to wide chars (utf16) for system call
        std::wstring wglobsrch(UTF_8::convert_to_utf16(globsrch));
        HANDLE hFile = FindFirstFileW(wglobsrch.c_str(), &data);
        // check if system called returned a result
        // ToDo: maybe we should print a debug message
        if (hFile == INVALID_HANDLE_VALUE) return -1;

        // read directory
        while (true)
        {
          try
          {
            // the system will report the filenames with wide chars (utf16)
            std::string entry = UTF_8::convert_from_utf16(data.cFileName);
            // check if file ending matches exactly
            if (!ends_with(entry, ".dll")) continue;
            // load the plugin and increase counter
            if (load_plugin(path + entry)) ++ loaded;
            // check if there should be more entries
            if (GetLastError() == ERROR_NO_MORE_FILES) break;
            // load next entry (check for return type)
            if (!FindNextFileW(hFile, &data)) break;
          }
          catch (...)
          {
            // report the error to the console (should not happen)
            // seems like we got strange data from the system call?
            std::cerr << "filename in plugin path has invalid utf8?" << std::endl;
          }
        }
      }
      catch (utf8::invalid_utf8)
      {
        // report the error to the console (should not happen)
        // implementors should make sure to provide valid utf8
        std::cerr << "plugin path contains invalid utf8" << std::endl;
      }

    #else

      DIR *dp;
      struct dirent *dirp;
      if((dp  = opendir(path.c_str())) == NULL) return -1;
      while ((dirp = readdir(dp)) != NULL) {
        #if __APPLE__
          if (!ends_with(dirp->d_name, ".dylib")) continue;
        #else
          if (!ends_with(dirp->d_name, ".so")) continue;
        #endif
        if (load_plugin(path + dirp->d_name)) ++ loaded;
      }
      closedir(dp);

    #endif
    return loaded;

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  void str_rtrim(std::string& str, const std::string& delimiters = " \f\n\r\t\v")
  {
    str.erase( str.find_last_not_of( delimiters ) + 1 );
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  PreValue::PreValue(ParserState pstate, bool d, bool e, bool i, Type ct)
  : Expression(pstate, d, e, i, ct)
  { }
  PreValue::PreValue(const PreValue* ptr)
  : Expression(ptr)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Value::Value(ParserState pstate, bool d, bool e, bool i, Type ct)
  : PreValue(pstate, d, e, i, ct)
  { }
  Value::Value(const Value* ptr)
  : PreValue(ptr)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  List::List(ParserState pstate, size_t size, enum Sass_Separator sep, bool argl, bool bracket)
  : Value(pstate),
    Vectorized<Expression_Obj>(size),
    separator_(sep),
    is_arglist_(argl),
    is_bracketed_(bracket),
    from_selector_(false)
  { concrete_type(LIST); }

  List::List(const List* ptr)
  : Value(ptr),
    Vectorized<Expression_Obj>(*ptr),
    separator_(ptr->separator_),
    is_arglist_(ptr->is_arglist_),
    is_bracketed_(ptr->is_bracketed_),
    from_selector_(ptr->from_selector_)
  { concrete_type(LIST); }

  size_t List::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<std::string>()(sep_string());
      hash_combine(hash_, std::hash<bool>()(is_bracketed()));
      for (size_t i = 0, L = length(); i < L; ++i)
        hash_combine(hash_, (elements()[i])->hash());
    }
    return hash_;
  }

  void List::set_delayed(bool delayed)
  {
    is_delayed(delayed);
    // don't set children
  }

  bool List::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<List>(&rhs)) {
      if (length() != r->length()) return false;
      if (separator() != r->separator()) return false;
      if (is_bracketed() != r->is_bracketed()) return false;
      for (size_t i = 0, L = length(); i < L; ++i) {
        auto rv = r->at(i);
        auto lv = this->at(i);
        if (!lv && rv) return false;
        else if (!rv && lv) return false;
        else if (*lv != *rv) return false;
      }
      return true;
    }
    return false;
  }

  size_t List::size() const {
    if (!is_arglist_) return length();
    // arglist expects a list of arguments
    // so we need to break before keywords
    for (size_t i = 0, L = length(); i < L; ++i) {
      Expression_Obj obj = this->at(i);
      if (Argument_Ptr arg = Cast<Argument>(obj)) {
        if (!arg->name().empty()) return i;
      }
    }
    return length();
  }


  Expression_Obj List::value_at_index(size_t i) {
    Expression_Obj obj = this->at(i);
    if (is_arglist_) {
      if (Argument_Ptr arg = Cast<Argument>(obj)) {
        return arg->value();
      } else {
        return obj;
      }
    } else {
      return obj;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Map::Map(ParserState pstate, size_t size)
  : Value(pstate),
    Hashed(size)
  { concrete_type(MAP); }

  Map::Map(const Map* ptr)
  : Value(ptr),
    Hashed(*ptr)
  { concrete_type(MAP); }

  bool Map::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Map>(&rhs)) {
      if (length() != r->length()) return false;
      for (auto key : keys()) {
        auto rv = r->at(key);
        auto lv = this->at(key);
        if (!lv && rv) return false;
        else if (!rv && lv) return false;
        else if (*lv != *rv) return false;
      }
      return true;
    }
    return false;
  }

  List_Obj Map::to_list(ParserState& pstate) {
    List_Obj ret = SASS_MEMORY_NEW(List, pstate, length(), SASS_COMMA);

    for (auto key : keys()) {
      List_Obj l = SASS_MEMORY_NEW(List, pstate, 2);
      l->append(key);
      l->append(at(key));
      ret->append(l);
    }

    return ret;
  }

  size_t Map::hash() const
  {
    if (hash_ == 0) {
      for (auto key : keys()) {
        hash_combine(hash_, key->hash());
        hash_combine(hash_, at(key)->hash());
      }
    }

    return hash_;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Binary_Expression::Binary_Expression(ParserState pstate,
                    Operand op, Expression_Obj lhs, Expression_Obj rhs)
  : PreValue(pstate), op_(op), left_(lhs), right_(rhs), hash_(0)
  { }

  Binary_Expression::Binary_Expression(const Binary_Expression* ptr)
  : PreValue(ptr),
    op_(ptr->op_),
    left_(ptr->left_),
    right_(ptr->right_),
    hash_(ptr->hash_)
  { }

  bool Binary_Expression::is_left_interpolant(void) const
  {
    return is_interpolant() || (left() && left()->is_left_interpolant());
  }
  bool Binary_Expression::is_right_interpolant(void) const
  {
    return is_interpolant() || (right() && right()->is_right_interpolant());
  }

  const std::string Binary_Expression::type_name()
  {
    return sass_op_to_name(optype());
  }

  const std::string Binary_Expression::separator()
  {
    return sass_op_separator(optype());
  }

  bool Binary_Expression::has_interpolant() const
  {
    return is_left_interpolant() ||
            is_right_interpolant();
  }

  void Binary_Expression::set_delayed(bool delayed)
  {
    right()->set_delayed(delayed);
    left()->set_delayed(delayed);
    is_delayed(delayed);
  }

  bool Binary_Expression::operator==(const Expression& rhs) const
  {
    if (auto m = Cast<Binary_Expression>(&rhs)) {
      return type() == m->type() &&
             *left() == *m->left() &&
             *right() == *m->right();
    }
    return false;
  }

  size_t Binary_Expression::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<size_t>()(optype());
      hash_combine(hash_, left()->hash());
      hash_combine(hash_, right()->hash());
    }
    return hash_;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Function::Function(ParserState pstate, Definition_Obj def, bool css)
  : Value(pstate), definition_(def), is_css_(css)
  { concrete_type(FUNCTION_VAL); }

  Function::Function(const Function* ptr)
  : Value(ptr), definition_(ptr->definition_), is_css_(ptr->is_css_)
  { concrete_type(FUNCTION_VAL); }

  bool Function::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Function>(&rhs)) {
      auto d1 = Cast<Definition>(definition());
      auto d2 = Cast<Definition>(r->definition());
      return d1 && d2 && d1 == d2 && is_css() == r->is_css();
    }
    return false;
  }

  std::string Function::name() {
    if (definition_) {
      return definition_->name();
    }
    return "";
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Function_Call::Function_Call(ParserState pstate, String_Obj n, Arguments_Obj args, void* cookie)
  : PreValue(pstate), sname_(n), arguments_(args), func_(), via_call_(false), cookie_(cookie), hash_(0)
  { concrete_type(FUNCTION); }
  Function_Call::Function_Call(ParserState pstate, String_Obj n, Arguments_Obj args, Function_Obj func)
  : PreValue(pstate), sname_(n), arguments_(args), func_(func), via_call_(false), cookie_(0), hash_(0)
  { concrete_type(FUNCTION); }
  Function_Call::Function_Call(ParserState pstate, String_Obj n, Arguments_Obj args)
  : PreValue(pstate), sname_(n), arguments_(args), via_call_(false), cookie_(0), hash_(0)
  { concrete_type(FUNCTION); }

  Function_Call::Function_Call(ParserState pstate, std::string n, Arguments_Obj args, void* cookie)
  : PreValue(pstate), sname_(SASS_MEMORY_NEW(String_Constant, pstate, n)), arguments_(args), func_(), via_call_(false), cookie_(cookie), hash_(0)
  { concrete_type(FUNCTION); }
  Function_Call::Function_Call(ParserState pstate, std::string n, Arguments_Obj args, Function_Obj func)
  : PreValue(pstate), sname_(SASS_MEMORY_NEW(String_Constant, pstate, n)), arguments_(args), func_(func), via_call_(false), cookie_(0), hash_(0)
  { concrete_type(FUNCTION); }
  Function_Call::Function_Call(ParserState pstate, std::string n, Arguments_Obj args)
  : PreValue(pstate), sname_(SASS_MEMORY_NEW(String_Constant, pstate, n)), arguments_(args), via_call_(false), cookie_(0), hash_(0)
  { concrete_type(FUNCTION); }

  Function_Call::Function_Call(const Function_Call* ptr)
  : PreValue(ptr),
    sname_(ptr->sname_),
    arguments_(ptr->arguments_),
    func_(ptr->func_),
    via_call_(ptr->via_call_),
    cookie_(ptr->cookie_),
    hash_(ptr->hash_)
  { concrete_type(FUNCTION); }

  bool Function_Call::operator==(const Expression& rhs) const
  {
    if (auto m = Cast<Function_Call>(&rhs)) {
      if (*sname() != *m->sname()) return false;
      if (arguments()->length() != m->arguments()->length()) return false;
      for (size_t i = 0, L = arguments()->length(); i < L; ++i)
        if (*arguments()->get(i) != *m->arguments()->get(i)) return false;
      return true;
    }
    return false;
  }

  size_t Function_Call::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<std::string>()(name());
      for (auto argument : arguments()->elements())
        hash_combine(hash_, argument->hash());
    }
    return hash_;
  }

  std::string Function_Call::name() const
  {
    return sname();
  }

  bool Function_Call::is_css() {
    if (func_) return func_->is_css();
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Variable::Variable(ParserState pstate, std::string n)
  : PreValue(pstate), name_(n)
  { concrete_type(VARIABLE); }

  Variable::Variable(const Variable* ptr)
  : PreValue(ptr), name_(ptr->name_)
  { concrete_type(VARIABLE); }

  bool Variable::operator==(const Expression& rhs) const
  {
    if (auto e = Cast<Variable>(&rhs)) {
      return name() == e->name();
    }
    return false;
  }

  size_t Variable::hash() const
  {
    return std::hash<std::string>()(name());
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Number::Number(ParserState pstate, double val, std::string u, bool zero)
  : Value(pstate),
    Units(),
    value_(val),
    zero_(zero),
    hash_(0)
  {
    size_t l = 0;
    size_t r;
    if (!u.empty()) {
      bool nominator = true;
      while (true) {
        r = u.find_first_of("*/", l);
        std::string unit(u.substr(l, r == std::string::npos ? r : r - l));
        if (!unit.empty()) {
          if (nominator) numerators.push_back(unit);
          else denominators.push_back(unit);
        }
        if (r == std::string::npos) break;
        // ToDo: should error for multiple slashes
        // if (!nominator && u[r] == '/') error(...)
        if (u[r] == '/')
          nominator = false;
        // strange math parsing?
        // else if (u[r] == '*')
        //  nominator = true;
        l = r + 1;
      }
    }
    concrete_type(NUMBER);
  }

  Number::Number(const Number* ptr)
  : Value(ptr),
    Units(ptr),
    value_(ptr->value_), zero_(ptr->zero_),
    hash_(ptr->hash_)
  { concrete_type(NUMBER); }

  // cancel out unnecessary units
  void Number::reduce()
  {
    // apply conversion factor
    value_ *= this->Units::reduce();
  }

  void Number::normalize()
  {
    // apply conversion factor
    value_ *= this->Units::normalize();
  }

  size_t Number::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<double>()(value_);
      for (const auto numerator : numerators)
        hash_combine(hash_, std::hash<std::string>()(numerator));
      for (const auto denominator : denominators)
        hash_combine(hash_, std::hash<std::string>()(denominator));
    }
    return hash_;
  }

  bool Number::operator== (const Expression& rhs) const
  {
    if (auto n = Cast<Number>(&rhs)) {
      return *this == *n;
    }
    return false;
  }

  bool Number::operator== (const Number& rhs) const
  {
    // unitless or only having one unit are equivalent (3.4)
    // therefore we need to reduce the units beforehand
    Number l(*this), r(rhs); l.reduce(); r.reduce();
    size_t lhs_units = l.numerators.size() + l.denominators.size();
    size_t rhs_units = r.numerators.size() + r.denominators.size();
    if (!lhs_units || !rhs_units) {
      return NEAR_EQUAL(l.value(), r.value());
    }
    // ensure both have same units
    l.normalize(); r.normalize();
    Units &lhs_unit = l, &rhs_unit = r;
    return lhs_unit == rhs_unit &&
      NEAR_EQUAL(l.value(), r.value());
  }

  bool Number::operator< (const Number& rhs) const
  {
    // unitless or only having one unit are equivalent (3.4)
    // therefore we need to reduce the units beforehand
    Number l(*this), r(rhs); l.reduce(); r.reduce();
    size_t lhs_units = l.numerators.size() + l.denominators.size();
    size_t rhs_units = r.numerators.size() + r.denominators.size();
    if (!lhs_units || !rhs_units) {
      return l.value() < r.value();
    }
    // ensure both have same units
    l.normalize(); r.normalize();
    Units &lhs_unit = l, &rhs_unit = r;
    if (!(lhs_unit == rhs_unit)) {
      /* ToDo: do we always get usefull backtraces? */
      throw Exception::IncompatibleUnits(rhs, *this);
    }
    if (lhs_unit == rhs_unit) {
      return l.value() < r.value();
    } else {
      return lhs_unit < rhs_unit;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Color::Color(ParserState pstate, double a, const std::string disp)
  : Value(pstate),
    disp_(disp), a_(a),
    hash_(0)
  { concrete_type(COLOR); }

  Color::Color(const Color* ptr)
  : Value(ptr->pstate()),
    // reset on copy
    disp_(""),
    a_(ptr->a_),
    hash_(ptr->hash_)
  { concrete_type(COLOR); }

  bool Color::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Color_RGBA>(&rhs)) {
      return *this == *r;
    }
    else if (auto r = Cast<Color_HSLA>(&rhs)) {
      return *this == *r;
    }
    else if (auto r = Cast<Color>(&rhs)) {
      return a_ == r->a();
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Color_RGBA::Color_RGBA(ParserState pstate, double r, double g, double b, double a, const std::string disp)
  : Color(pstate, a, disp),
    r_(r), g_(g), b_(b)
  { concrete_type(COLOR); }

  Color_RGBA::Color_RGBA(const Color_RGBA* ptr)
  : Color(ptr),
    r_(ptr->r_),
    g_(ptr->g_),
    b_(ptr->b_)
  { concrete_type(COLOR); }

  bool Color_RGBA::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Color_RGBA>(&rhs)) {
      return r_ == r->r() &&
             g_ == r->g() &&
             b_ == r->b() &&
             a_ == r->a();
    }
    return false;
  }

  size_t Color_RGBA::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<std::string>()("RGBA");
      hash_combine(hash_, std::hash<double>()(a_));
      hash_combine(hash_, std::hash<double>()(r_));
      hash_combine(hash_, std::hash<double>()(g_));
      hash_combine(hash_, std::hash<double>()(b_));
    }
    return hash_;
  }

  Color_HSLA_Ptr Color_RGBA::toHSLA(bool copy)
  {

    // Algorithm from http://en.wikipedia.org/wiki/wHSL_and_HSV#Conversion_from_RGB_to_HSL_or_HSV
    double r = r_ / 255.0;
    double g = g_ / 255.0;
    double b = b_ / 255.0;

    double max = std::max(r, std::max(g, b));
    double min = std::min(r, std::min(g, b));
    double delta = max - min;

    double h = 0;
    double s;
    double l = (max + min) / 2.0;

    if (NEAR_EQUAL(max, min)) {
      h = s = 0; // achromatic
    }
    else {
      if (l < 0.5) s = delta / (max + min);
      else         s = delta / (2.0 - max - min);

      if      (r == max) h = (g - b) / delta + (g < b ? 6 : 0);
      else if (g == max) h = (b - r) / delta + 2;
      else if (b == max) h = (r - g) / delta + 4;
    }

    // HSL hsl_struct;
    h = h * 60;
    s = s * 100;
    l = l * 100;

    return SASS_MEMORY_NEW(Color_HSLA,
      pstate(), h, s, l, a(), ""
    );
  }

  Color_RGBA_Ptr Color_RGBA::toRGBA(bool copy)
  {
    return copy ? SASS_MEMORY_COPY(this) : this;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Color_HSLA::Color_HSLA(ParserState pstate, double h, double s, double l, double a, const std::string disp)
  : Color(pstate, a, disp),
    h_(absmod(h, 360.0)),
    s_(clip(s, 0.0, 100.0)),
    l_(clip(l, 0.0, 100.0))
    // hash_(0)
  { concrete_type(COLOR); }

  Color_HSLA::Color_HSLA(const Color_HSLA* ptr)
  : Color(ptr),
    h_(ptr->h_),
    s_(ptr->s_),
    l_(ptr->l_)
    // hash_(ptr->hash_)
  { concrete_type(COLOR); }

  bool Color_HSLA::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Color_HSLA>(&rhs)) {
      return h_ == r->h() &&
             s_ == r->s() &&
             l_ == r->l() &&
             a_ == r->a();
    }
    return false;
  }

  size_t Color_HSLA::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<std::string>()("HSLA");
      hash_combine(hash_, std::hash<double>()(a_));
      hash_combine(hash_, std::hash<double>()(h_));
      hash_combine(hash_, std::hash<double>()(s_));
      hash_combine(hash_, std::hash<double>()(l_));
    }
    return hash_;
  }

  // hue to RGB helper function
  double h_to_rgb(double m1, double m2, double h)
  {
    h = absmod(h, 1.0);
    if (h*6.0 < 1) return m1 + (m2 - m1)*h*6;
    if (h*2.0 < 1) return m2;
    if (h*3.0 < 2) return m1 + (m2 - m1) * (2.0/3.0 - h)*6;
    return m1;
  }

  Color_RGBA_Ptr Color_HSLA::toRGBA(bool copy)
  {

    double h = absmod(h_ / 360.0, 1.0);
    double s = clip(s_ / 100.0, 0.0, 1.0);
    double l = clip(l_ / 100.0, 0.0, 1.0);

    // Algorithm from the CSS3 spec: http://www.w3.org/TR/css3-color/#hsl-color.
    double m2;
    if (l <= 0.5) m2 = l*(s+1.0);
    else m2 = (l+s)-(l*s);
    double m1 = (l*2.0)-m2;
    // round the results -- consider moving this into the Color constructor
    double r = (h_to_rgb(m1, m2, h + 1.0/3.0) * 255.0);
    double g = (h_to_rgb(m1, m2, h) * 255.0);
    double b = (h_to_rgb(m1, m2, h - 1.0/3.0) * 255.0);

    return SASS_MEMORY_NEW(Color_RGBA,
      pstate(), r, g, b, a(), ""
    );
  }

  Color_HSLA_Ptr Color_HSLA::toHSLA(bool copy)
  {
    return copy ? SASS_MEMORY_COPY(this) : this;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Custom_Error::Custom_Error(ParserState pstate, std::string msg)
  : Value(pstate), message_(msg)
  { concrete_type(C_ERROR); }

  Custom_Error::Custom_Error(const Custom_Error* ptr)
  : Value(ptr), message_(ptr->message_)
  { concrete_type(C_ERROR); }

  bool Custom_Error::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Custom_Error>(&rhs)) {
      return message() == r->message();
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Custom_Warning::Custom_Warning(ParserState pstate, std::string msg)
  : Value(pstate), message_(msg)
  { concrete_type(C_WARNING); }

  Custom_Warning::Custom_Warning(const Custom_Warning* ptr)
  : Value(ptr), message_(ptr->message_)
  { concrete_type(C_WARNING); }

  bool Custom_Warning::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Custom_Warning>(&rhs)) {
      return message() == r->message();
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Boolean::Boolean(ParserState pstate, bool val)
  : Value(pstate), value_(val),
    hash_(0)
  { concrete_type(BOOLEAN); }

  Boolean::Boolean(const Boolean* ptr)
  : Value(ptr),
    value_(ptr->value_),
    hash_(ptr->hash_)
  { concrete_type(BOOLEAN); }

 bool Boolean::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<Boolean>(&rhs)) {
      return (value() == r->value());
    }
    return false;
  }

  size_t Boolean::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<bool>()(value_);
    }
    return hash_;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  String::String(ParserState pstate, bool delayed)
  : Value(pstate, delayed)
  { concrete_type(STRING); }
  String::String(const String* ptr)
  : Value(ptr)
  { concrete_type(STRING); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  String_Schema::String_Schema(ParserState pstate, size_t size, bool css)
  : String(pstate), Vectorized<PreValue_Obj>(size), css_(css), hash_(0)
  { concrete_type(STRING); }

  String_Schema::String_Schema(const String_Schema* ptr)
  : String(ptr),
    Vectorized<PreValue_Obj>(*ptr),
    css_(ptr->css_),
    hash_(ptr->hash_)
  { concrete_type(STRING); }

  void String_Schema::rtrim()
  {
    if (!empty()) {
      if (String_Ptr str = Cast<String>(last())) str->rtrim();
    }
  }

  bool String_Schema::is_left_interpolant(void) const
  {
    return length() && first()->is_left_interpolant();
  }
  bool String_Schema::is_right_interpolant(void) const
  {
    return length() && last()->is_right_interpolant();
  }

  bool String_Schema::operator== (const Expression& rhs) const
  {
    if (auto r = Cast<String_Schema>(&rhs)) {
      if (length() != r->length()) return false;
      for (size_t i = 0, L = length(); i < L; ++i) {
        auto rv = (*r)[i];
        auto lv = (*this)[i];
        if (*lv != *rv) return false;
      }
      return true;
    }
    return false;
  }

  bool String_Schema::has_interpolants()
  {
    for (auto el : elements()) {
      if (el->is_interpolant()) return true;
    }
    return false;
  }

  size_t String_Schema::hash() const
  {
    if (hash_ == 0) {
      for (auto string : elements())
        hash_combine(hash_, string->hash());
    }
    return hash_;
  }

  void String_Schema::set_delayed(bool delayed)
  {
    is_delayed(delayed);
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  String_Constant::String_Constant(ParserState pstate, std::string val, bool css)
  : String(pstate), quote_mark_(0), can_compress_whitespace_(false), value_(read_css_string(val, css)), hash_(0)
  { }
  String_Constant::String_Constant(ParserState pstate, const char* beg, bool css)
  : String(pstate), quote_mark_(0), can_compress_whitespace_(false), value_(read_css_string(std::string(beg), css)), hash_(0)
  { }
  String_Constant::String_Constant(ParserState pstate, const char* beg, const char* end, bool css)
  : String(pstate), quote_mark_(0), can_compress_whitespace_(false), value_(read_css_string(std::string(beg, end-beg), css)), hash_(0)
  { }
  String_Constant::String_Constant(ParserState pstate, const Token& tok, bool css)
  : String(pstate), quote_mark_(0), can_compress_whitespace_(false), value_(read_css_string(std::string(tok.begin, tok.end), css)), hash_(0)
  { }

  String_Constant::String_Constant(const String_Constant* ptr)
  : String(ptr),
    quote_mark_(ptr->quote_mark_),
    can_compress_whitespace_(ptr->can_compress_whitespace_),
    value_(ptr->value_),
    hash_(ptr->hash_)
  { }

  bool String_Constant::is_invisible() const {
    return value_.empty() && quote_mark_ == 0;
  }

  bool String_Constant::operator== (const Expression& rhs) const
  {
    if (auto qstr = Cast<String_Quoted>(&rhs)) {
      return value() == qstr->value();
    } else if (auto cstr = Cast<String_Constant>(&rhs)) {
      return value() == cstr->value();
    }
    return false;
  }

  std::string String_Constant::inspect() const
  {
    return quote(value_, '*');
  }

  void String_Constant::rtrim()
  {
    str_rtrim(value_);
  }

  size_t String_Constant::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<std::string>()(value_);
    }
    return hash_;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  String_Quoted::String_Quoted(ParserState pstate, std::string val, char q,
    bool keep_utf8_escapes, bool skip_unquoting,
    bool strict_unquoting, bool css)
  : String_Constant(pstate, val, css)
  {
    if (skip_unquoting == false) {
      value_ = unquote(value_, &quote_mark_, keep_utf8_escapes, strict_unquoting);
    }
    if (q && quote_mark_) quote_mark_ = q;
  }

  String_Quoted::String_Quoted(const String_Quoted* ptr)
  : String_Constant(ptr)
  { }

  bool String_Quoted::operator== (const Expression& rhs) const
  {
    if (auto qstr = Cast<String_Quoted>(&rhs)) {
      return value() == qstr->value();
    } else if (auto cstr = Cast<String_Constant>(&rhs)) {
      return value() == cstr->value();
    }
    return false;
  }

  std::string String_Quoted::inspect() const
  {
    return quote(value_, '*');
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Null::Null(ParserState pstate)
  : Value(pstate)
  { concrete_type(NULL_VAL); }

  Null::Null(const Null* ptr) : Value(ptr)
  { concrete_type(NULL_VAL); }

  bool Null::operator== (const Expression& rhs) const
  {
    return Cast<Null>(&rhs) != NULL;
  }

  size_t Null::hash() const
  {
    return -1;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Parent_Reference::Parent_Reference(ParserState pstate)
  : Value(pstate)
  { concrete_type(PARENT); }

  Parent_Reference::Parent_Reference(const Parent_Reference* ptr)
  : Value(ptr)
  { concrete_type(PARENT); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  IMPLEMENT_AST_OPERATORS(List);
  IMPLEMENT_AST_OPERATORS(Map);
  IMPLEMENT_AST_OPERATORS(Binary_Expression);
  IMPLEMENT_AST_OPERATORS(Function);
  IMPLEMENT_AST_OPERATORS(Function_Call);
  IMPLEMENT_AST_OPERATORS(Variable);
  IMPLEMENT_AST_OPERATORS(Number);
  IMPLEMENT_AST_OPERATORS(Color_RGBA);
  IMPLEMENT_AST_OPERATORS(Color_HSLA);
  IMPLEMENT_AST_OPERATORS(Custom_Error);
  IMPLEMENT_AST_OPERATORS(Custom_Warning);
  IMPLEMENT_AST_OPERATORS(Boolean);
  IMPLEMENT_AST_OPERATORS(String_Schema);
  IMPLEMENT_AST_OPERATORS(String_Constant);
  IMPLEMENT_AST_OPERATORS(String_Quoted);
  IMPLEMENT_AST_OPERATORS(Null);
  IMPLEMENT_AST_OPERATORS(Parent_Reference);

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

}

namespace Sass {

  namespace Functions {

    // features
    static std::set<std::string> features {
      "global-variable-shadowing",
      "extend-selector-pseudoclass",
      "at-error",
      "units-level-3",
      "custom-property"
    };

    //////////////////////////
    // INTROSPECTION FUNCTIONS
    //////////////////////////

    Signature type_of_sig = "type-of($value)";
    BUILT_IN(type_of)
    {
      Expression_Ptr v = ARG("$value", Expression);
      return SASS_MEMORY_NEW(String_Quoted, pstate, v->type());
    }

    Signature variable_exists_sig = "variable-exists($name)";
    BUILT_IN(variable_exists)
    {
      std::string s = Util::normalize_underscores(unquote(ARG("$name", String_Constant)->value()));

      if(d_env.has("$"+s)) {
        return SASS_MEMORY_NEW(Boolean, pstate, true);
      }
      else {
        return SASS_MEMORY_NEW(Boolean, pstate, false);
      }
    }

    Signature global_variable_exists_sig = "global-variable-exists($name)";
    BUILT_IN(global_variable_exists)
    {
      std::string s = Util::normalize_underscores(unquote(ARG("$name", String_Constant)->value()));

      if(d_env.has_global("$"+s)) {
        return SASS_MEMORY_NEW(Boolean, pstate, true);
      }
      else {
        return SASS_MEMORY_NEW(Boolean, pstate, false);
      }
    }

    Signature function_exists_sig = "function-exists($name)";
    BUILT_IN(function_exists)
    {
      String_Constant_Ptr ss = Cast<String_Constant>(env["$name"]);
      if (!ss) {
        error("$name: " + (env["$name"]->to_string()) + " is not a string for `function-exists'", pstate, traces);
      }

      std::string name = Util::normalize_underscores(unquote(ss->value()));

      if(d_env.has_global(name+"[f]")) {
        return SASS_MEMORY_NEW(Boolean, pstate, true);
      }
      else {
        return SASS_MEMORY_NEW(Boolean, pstate, false);
      }
    }

    Signature mixin_exists_sig = "mixin-exists($name)";
    BUILT_IN(mixin_exists)
    {
      std::string s = Util::normalize_underscores(unquote(ARG("$name", String_Constant)->value()));

      if(d_env.has_global(s+"[m]")) {
        return SASS_MEMORY_NEW(Boolean, pstate, true);
      }
      else {
        return SASS_MEMORY_NEW(Boolean, pstate, false);
      }
    }

    Signature feature_exists_sig = "feature-exists($name)";
    BUILT_IN(feature_exists)
    {
      std::string s = unquote(ARG("$name", String_Constant)->value());

      if(features.find(s) == features.end()) {
        return SASS_MEMORY_NEW(Boolean, pstate, false);
      }
      else {
        return SASS_MEMORY_NEW(Boolean, pstate, true);
      }
    }

    Signature call_sig = "call($name, $args...)";
    BUILT_IN(call)
    {
      std::string name;
      Function_Ptr ff = Cast<Function>(env["$name"]);
      String_Constant_Ptr ss = Cast<String_Constant>(env["$name"]);

      if (ss) {
        name = Util::normalize_underscores(unquote(ss->value()));
        std::cerr << "DEPRECATION WARNING: ";
        std::cerr << "Passing a string to call() is deprecated and will be illegal" << std::endl;
        std::cerr << "in Sass 4.0. Use call(get-function(" + quote(name) + ")) instead." << std::endl;
        std::cerr << std::endl;
      } else if (ff) {
        name = ff->name();
      }

      List_Obj arglist = SASS_MEMORY_COPY(ARG("$args", List));

      Arguments_Obj args = SASS_MEMORY_NEW(Arguments, pstate);
      // std::string full_name(name + "[f]");
      // Definition_Ptr def = d_env.has(full_name) ? Cast<Definition>((d_env)[full_name]) : 0;
      // Parameters_Ptr params = def ? def->parameters() : 0;
      // size_t param_size = params ? params->length() : 0;
      for (size_t i = 0, L = arglist->length(); i < L; ++i) {
        Expression_Obj expr = arglist->value_at_index(i);
        // if (params && params->has_rest_parameter()) {
        //   Parameter_Obj p = param_size > i ? (*params)[i] : 0;
        //   List_Ptr list = Cast<List>(expr);
        //   if (list && p && !p->is_rest_parameter()) expr = (*list)[0];
        // }
        if (arglist->is_arglist()) {
          Expression_Obj obj = arglist->at(i);
          Argument_Obj arg = (Argument_Ptr) obj.ptr(); // XXX
          args->append(SASS_MEMORY_NEW(Argument,
                                       pstate,
                                       expr,
                                       arg ? arg->name() : "",
                                       arg ? arg->is_rest_argument() : false,
                                       arg ? arg->is_keyword_argument() : false));
        } else {
          args->append(SASS_MEMORY_NEW(Argument, pstate, expr));
        }
      }
      Function_Call_Obj func = SASS_MEMORY_NEW(Function_Call, pstate, name, args);
      Expand expand(ctx, &d_env, &selector_stack);
      func->via_call(true); // calc invoke is allowed
      if (ff) func->func(ff);
      return Cast<PreValue>(func->perform(&expand.eval));
    }

    ////////////////////
    // BOOLEAN FUNCTIONS
    ////////////////////

    Signature not_sig = "not($value)";
    BUILT_IN(sass_not)
    {
      return SASS_MEMORY_NEW(Boolean, pstate, ARG("$value", Expression)->is_false());
    }

    Signature if_sig = "if($condition, $if-true, $if-false)";
    // BUILT_IN(sass_if)
    // { return ARG("$condition", Expression)->is_false() ? ARG("$if-false", Expression) : ARG("$if-true", Expression); }
    BUILT_IN(sass_if)
    {
      Expand expand(ctx, &d_env, &selector_stack);
      Expression_Obj cond = ARG("$condition", Expression)->perform(&expand.eval);
      bool is_true = !cond->is_false();
      Expression_Obj res = ARG(is_true ? "$if-true" : "$if-false", Expression);
      Value_Obj qwe = Cast<Value>(res->perform(&expand.eval));
      // res = res->perform(&expand.eval.val_eval);
      qwe->set_delayed(false); // clone?
      return qwe.detach();
    }

    //////////////////////////
    // MISCELLANEOUS FUNCTIONS
    //////////////////////////

    // value.check_deprecated_interp if value.is_a?(Sass::Script::Value::String)
    // unquoted_string(value.to_sass)

    Signature inspect_sig = "inspect($value)";
    BUILT_IN(inspect)
    {
      Expression_Ptr v = ARG("$value", Expression);
      if (v->concrete_type() == Expression::NULL_VAL) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "null");
      } else if (v->concrete_type() == Expression::BOOLEAN && v->is_false()) {
        return SASS_MEMORY_NEW(String_Quoted, pstate, "false");
      } else if (v->concrete_type() == Expression::STRING) {
        return Cast<String>(v);
      } else {
        // ToDo: fix to_sass for nested parentheses
        Sass_Output_Style old_style;
        old_style = ctx.c_options.output_style;
        ctx.c_options.output_style = TO_SASS;
        Emitter emitter(ctx.c_options);
        Inspect i(emitter);
        i.in_declaration = false;
        v->perform(&i);
        ctx.c_options.output_style = old_style;
        return SASS_MEMORY_NEW(String_Quoted, pstate, i.get_buffer());
      }
      // return v;
    }

    Signature content_exists_sig = "content-exists()";
    BUILT_IN(content_exists)
    {
      if (!d_env.has_global("is_in_mixin")) {
        error("Cannot call content-exists() except within a mixin.", pstate, traces);
      }
      return SASS_MEMORY_NEW(Boolean, pstate, d_env.has_lexical("@content[m]"));
    }

    Signature get_function_sig = "get-function($name, $css: false)";
    BUILT_IN(get_function)
    {
      String_Constant_Ptr ss = Cast<String_Constant>(env["$name"]);
      if (!ss) {
        error("$name: " + (env["$name"]->to_string()) + " is not a string for `get-function'", pstate, traces);
      }

      std::string name = Util::normalize_underscores(unquote(ss->value()));
      std::string full_name = name + "[f]";

      Boolean_Obj css = ARG("$css", Boolean);
      if (!css->is_false()) {
        Definition_Ptr def = SASS_MEMORY_NEW(Definition,
                                         pstate,
                                         name,
                                         SASS_MEMORY_NEW(Parameters, pstate),
                                         SASS_MEMORY_NEW(Block, pstate, 0, false),
                                         Definition::FUNCTION);
        return SASS_MEMORY_NEW(Function, pstate, def, true);
      }


      if (!d_env.has_global(full_name)) {
        error("Function not found: " + name, pstate, traces);
      }

      Definition_Ptr def = Cast<Definition>(d_env[full_name]);
      return SASS_MEMORY_NEW(Function, pstate, def, false);
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  Listize::Listize()
  {  }

  Expression_Ptr Listize::operator()(Selector_List_Ptr sel)
  {
    List_Obj l = SASS_MEMORY_NEW(List, sel->pstate(), sel->length(), SASS_COMMA);
    l->from_selector(true);
    for (size_t i = 0, L = sel->length(); i < L; ++i) {
      if (!sel->at(i)) continue;
      l->append(sel->at(i)->perform(this));
    }
    if (l->length()) return l.detach();
    return SASS_MEMORY_NEW(Null, l->pstate());
  }

  Expression_Ptr Listize::operator()(Compound_Selector_Ptr sel)
  {
    std::string str;
    for (size_t i = 0, L = sel->length(); i < L; ++i) {
      Expression_Ptr e = (*sel)[i]->perform(this);
      if (e) str += e->to_string();
    }
    return SASS_MEMORY_NEW(String_Quoted, sel->pstate(), str);
  }

  Expression_Ptr Listize::operator()(Complex_Selector_Ptr sel)
  {
    List_Obj l = SASS_MEMORY_NEW(List, sel->pstate(), 2);
    l->from_selector(true);
    Compound_Selector_Obj head = sel->head();
    if (head && !head->is_empty_reference())
    {
      Expression_Ptr hh = head->perform(this);
      if (hh) l->append(hh);
    }

    std::string reference = ! sel->reference() ? ""
      : sel->reference()->to_string();
    switch(sel->combinator())
    {
      case Complex_Selector::PARENT_OF:
        l->append(SASS_MEMORY_NEW(String_Quoted, sel->pstate(), ">"));
      break;
      case Complex_Selector::ADJACENT_TO:
        l->append(SASS_MEMORY_NEW(String_Quoted, sel->pstate(), "+"));
      break;
      case Complex_Selector::REFERENCE:
        l->append(SASS_MEMORY_NEW(String_Quoted, sel->pstate(), "/" + reference + "/"));
      break;
      case Complex_Selector::PRECEDES:
        l->append(SASS_MEMORY_NEW(String_Quoted, sel->pstate(), "~"));
      break;
      case Complex_Selector::ANCESTOR_OF:
      break;
      default: break;
    }

    Complex_Selector_Obj tail = sel->tail();
    if (tail)
    {
      Expression_Obj tt = tail->perform(this);
      if (List_Ptr ls = Cast<List>(tt))
      { l->concat(ls); }
    }
    if (l->length() == 0) return 0;
    return l.detach();
  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  // convert value from C++ side to C-API
  union Sass_Value* ast_node_to_sass_value (const Expression_Ptr val)
  {
    if (val->concrete_type() == Expression::NUMBER)
    {
      Number_Ptr_Const res = Cast<Number>(val);
      return sass_make_number(res->value(), res->unit().c_str());
    }
    else if (val->concrete_type() == Expression::COLOR)
    {
      // ToDo: allow to also use HSLA colors!!
      Color_RGBA_Obj col = Cast<Color>(val)->toRGBA();
      return sass_make_color(col->r(), col->g(), col->b(), col->a());
    }
    else if (val->concrete_type() == Expression::LIST)
    {
      List_Ptr_Const l = Cast<List>(val);
      union Sass_Value* list = sass_make_list(l->size(), l->separator(), l->is_bracketed());
      for (size_t i = 0, L = l->length(); i < L; ++i) {
        Expression_Obj obj = l->at(i);
        auto val = ast_node_to_sass_value(obj);
        sass_list_set_value(list, i, val);
      }
      return list;
    }
    else if (val->concrete_type() == Expression::MAP)
    {
      Map_Ptr_Const m = Cast<Map>(val);
      union Sass_Value* map = sass_make_map(m->length());
      size_t i = 0; for (Expression_Obj key : m->keys()) {
        sass_map_set_key(map, i, ast_node_to_sass_value(key));
        sass_map_set_value(map, i, ast_node_to_sass_value(m->at(key)));
        ++ i;
      }
      return map;
    }
    else if (val->concrete_type() == Expression::NULL_VAL)
    {
      return sass_make_null();
    }
    else if (val->concrete_type() == Expression::BOOLEAN)
    {
      Boolean_Ptr_Const res = Cast<Boolean>(val);
      return sass_make_boolean(res->value());
    }
    else if (val->concrete_type() == Expression::STRING)
    {
      if (String_Quoted_Ptr_Const qstr = Cast<String_Quoted>(val))
      {
        return sass_make_qstring(qstr->value().c_str());
      }
      else if (String_Constant_Ptr_Const cstr = Cast<String_Constant>(val))
      {
        return sass_make_string(cstr->value().c_str());
      }
    }
    return sass_make_error("unknown sass value type");
  }

  // convert value from C-API to C++ side
  Value_Ptr sass_value_to_ast_node (const union Sass_Value* val)
  {
    switch (sass_value_get_tag(val)) {
      case SASS_NUMBER:
        return SASS_MEMORY_NEW(Number,
                               ParserState("[C-VALUE]"),
                               sass_number_get_value(val),
                               sass_number_get_unit(val));
      case SASS_BOOLEAN:
        return SASS_MEMORY_NEW(Boolean,
                               ParserState("[C-VALUE]"),
                               sass_boolean_get_value(val));
      case SASS_COLOR:
        // ToDo: allow to also use HSLA colors!!
        return SASS_MEMORY_NEW(Color_RGBA,
                               ParserState("[C-VALUE]"),
                               sass_color_get_r(val),
                               sass_color_get_g(val),
                               sass_color_get_b(val),
                               sass_color_get_a(val));
      case SASS_STRING:
        if (sass_string_is_quoted(val)) {
          return SASS_MEMORY_NEW(String_Quoted,
                                 ParserState("[C-VALUE]"),
                                 sass_string_get_value(val));
        }
        return SASS_MEMORY_NEW(String_Constant,
                                 ParserState("[C-VALUE]"),
                                 sass_string_get_value(val));
      case SASS_LIST: {
        List_Ptr l = SASS_MEMORY_NEW(List,
                                  ParserState("[C-VALUE]"),
                                  sass_list_get_length(val),
                                  sass_list_get_separator(val));
        for (size_t i = 0, L = sass_list_get_length(val); i < L; ++i) {
          l->append(sass_value_to_ast_node(sass_list_get_value(val, i)));
        }
        l->is_bracketed(sass_list_get_is_bracketed(val));
        return l;
      }
      case SASS_MAP: {
        Map_Ptr m = SASS_MEMORY_NEW(Map, ParserState("[C-VALUE]"));
        for (size_t i = 0, L = sass_map_get_length(val); i < L; ++i) {
          *m << std::make_pair(
            sass_value_to_ast_node(sass_map_get_key(val, i)),
            sass_value_to_ast_node(sass_map_get_value(val, i)));
        }
        return m;
      }
      case SASS_NULL:
        return SASS_MEMORY_NEW(Null, ParserState("[C-VALUE]"));
      case SASS_ERROR:
        return SASS_MEMORY_NEW(Custom_Error,
                               ParserState("[C-VALUE]"),
                               sass_error_get_message(val));
      case SASS_WARNING:
        return SASS_MEMORY_NEW(Custom_Warning,
                               ParserState("[C-VALUE]"),
                               sass_warning_get_message(val));
      default: break;
    }
    return 0;
  }

}

namespace Sass {

  static Null sass_null(ParserState("null"));

  const char* sass_op_to_name(enum Sass_OP op) {
    switch (op) {
      case AND: return "and";
      case OR: return "or";
      case EQ: return "eq";
      case NEQ: return "neq";
      case GT: return "gt";
      case GTE: return "gte";
      case LT: return "lt";
      case LTE: return "lte";
      case ADD: return "plus";
      case SUB: return "minus";
      case MUL: return "times";
      case DIV: return "div";
      case MOD: return "mod";
      // this is only used internally!
      case NUM_OPS: return "[OPS]";
      default: return "invalid";
    }
  }

  const char* sass_op_separator(enum Sass_OP op) {
    switch (op) {
      case AND: return "&&";
      case OR: return "||";
      case EQ: return "==";
      case NEQ: return "!=";
      case GT: return ">";
      case GTE: return ">=";
      case LT: return "<";
      case LTE: return "<=";
      case ADD: return "+";
      case SUB: return "-";
      case MUL: return "*";
      case DIV: return "/";
      case MOD: return "%";
      // this is only used internally!
      case NUM_OPS: return "[OPS]";
      default: return "invalid";
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  void AST_Node::update_pstate(const ParserState& pstate)
  {
    pstate_.offset += pstate - pstate_ + pstate.offset;
  }

  const std::string AST_Node::to_string(Sass_Inspect_Options opt) const
  {
    Sass_Output_Options out(opt);
    Emitter emitter(out);
    Inspect i(emitter);
    i.in_declaration = true;
    // ToDo: inspect should be const
    const_cast<AST_Node_Ptr>(this)->perform(&i);
    return i.get_buffer();
  }

  const std::string AST_Node::to_string() const
  {
    return to_string({ NESTED, 5 });
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Expression_Obj Hashed::at(Expression_Obj k) const
  {
    if (elements_.count(k))
    { return elements_.at(k); }
    else { return {}; }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Statement::Statement(ParserState pstate, Type st, size_t t)
  : AST_Node(pstate), statement_type_(st), tabs_(t), group_end_(false)
  { }
  Statement::Statement(const Statement* ptr)
  : AST_Node(ptr),
    statement_type_(ptr->statement_type_),
    tabs_(ptr->tabs_),
    group_end_(ptr->group_end_)
  { }

  bool Statement::bubbles()
  {
    return false;
  }

  bool Statement::has_content()
  {
    return statement_type_ == CONTENT;
  }

  bool Statement::is_invisible() const
  {
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Block::Block(ParserState pstate, size_t s, bool r)
  : Statement(pstate),
    Vectorized<Statement_Obj>(s),
    is_root_(r)
  { }
  Block::Block(const Block* ptr)
  : Statement(ptr),
    Vectorized<Statement_Obj>(*ptr),
    is_root_(ptr->is_root_)
  { }

  bool Block::has_content()
  {
    for (size_t i = 0, L = elements().size(); i < L; ++i) {
      if (elements()[i]->has_content()) return true;
    }
    return Statement::has_content();
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Has_Block::Has_Block(ParserState pstate, Block_Obj b)
  : Statement(pstate), block_(b)
  { }
  Has_Block::Has_Block(const Has_Block* ptr)
  : Statement(ptr), block_(ptr->block_)
  { }

  bool Has_Block::has_content()
  {
    return (block_ && block_->has_content()) || Statement::has_content();
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Ruleset::Ruleset(ParserState pstate, Selector_List_Obj s, Block_Obj b)
  : Has_Block(pstate, b), selector_(s), is_root_(false)
  { statement_type(RULESET); }
  Ruleset::Ruleset(const Ruleset* ptr)
  : Has_Block(ptr),
    selector_(ptr->selector_),
    is_root_(ptr->is_root_)
  { statement_type(RULESET); }

  bool Ruleset::is_invisible() const {
    if (Selector_List_Ptr sl = Cast<Selector_List>(selector())) {
      for (size_t i = 0, L = sl->length(); i < L; ++i)
        if (!(*sl)[i]->has_placeholder()) return false;
    }
    return true;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Bubble::Bubble(ParserState pstate, Statement_Obj n, Statement_Obj g, size_t t)
  : Statement(pstate, Statement::BUBBLE, t), node_(n), group_end_(g == 0)
  { }
  Bubble::Bubble(const Bubble* ptr)
  : Statement(ptr),
    node_(ptr->node_),
    group_end_(ptr->group_end_)
  { }

  bool Bubble::bubbles()
  {
    return true;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Trace::Trace(ParserState pstate, std::string n, Block_Obj b, char type)
  : Has_Block(pstate, b), type_(type), name_(n)
  { }
  Trace::Trace(const Trace* ptr)
  : Has_Block(ptr),
    type_(ptr->type_),
    name_(ptr->name_)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Media_Block::Media_Block(ParserState pstate, List_Obj mqs, Block_Obj b)
  : Has_Block(pstate, b), media_queries_(mqs)
  { statement_type(MEDIA); }
  Media_Block::Media_Block(const Media_Block* ptr)
  : Has_Block(ptr), media_queries_(ptr->media_queries_)
  { statement_type(MEDIA); }

  bool Media_Block::is_invisible() const {
    for (size_t i = 0, L = block()->length(); i < L; ++i) {
      Statement_Obj stm = block()->at(i);
      if (!stm->is_invisible()) return false;
    }
    return true;
  }

  bool Media_Block::bubbles()
  {
    return true;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Directive::Directive(ParserState pstate, std::string kwd, Selector_List_Obj sel, Block_Obj b, Expression_Obj val)
  : Has_Block(pstate, b), keyword_(kwd), selector_(sel), value_(val) // set value manually if needed
  { statement_type(DIRECTIVE); }
  Directive::Directive(const Directive* ptr)
  : Has_Block(ptr),
    keyword_(ptr->keyword_),
    selector_(ptr->selector_),
    value_(ptr->value_) // set value manually if needed
  { statement_type(DIRECTIVE); }

  bool Directive::bubbles() { return is_keyframes() || is_media(); }

  bool Directive::is_media() {
    return keyword_.compare("@-webkit-media") == 0 ||
            keyword_.compare("@-moz-media") == 0 ||
            keyword_.compare("@-o-media") == 0 ||
            keyword_.compare("@media") == 0;
  }
  bool Directive::is_keyframes() {
    return keyword_.compare("@-webkit-keyframes") == 0 ||
            keyword_.compare("@-moz-keyframes") == 0 ||
            keyword_.compare("@-o-keyframes") == 0 ||
            keyword_.compare("@keyframes") == 0;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Keyframe_Rule::Keyframe_Rule(ParserState pstate, Block_Obj b)
  : Has_Block(pstate, b), name_()
  { statement_type(KEYFRAMERULE); }
  Keyframe_Rule::Keyframe_Rule(const Keyframe_Rule* ptr)
  : Has_Block(ptr), name_(ptr->name_)
  { statement_type(KEYFRAMERULE); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Declaration::Declaration(ParserState pstate, String_Obj prop, Expression_Obj val, bool i, bool c, Block_Obj b)
  : Has_Block(pstate, b), property_(prop), value_(val), is_important_(i), is_custom_property_(c), is_indented_(false)
  { statement_type(DECLARATION); }
  Declaration::Declaration(const Declaration* ptr)
  : Has_Block(ptr),
    property_(ptr->property_),
    value_(ptr->value_),
    is_important_(ptr->is_important_),
    is_custom_property_(ptr->is_custom_property_),
    is_indented_(ptr->is_indented_)
  { statement_type(DECLARATION); }

  bool Declaration::is_invisible() const
  {
    if (is_custom_property()) return false;
    return !(value_ && !Cast<Null>(value_));
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Assignment::Assignment(ParserState pstate, std::string var, Expression_Obj val, bool is_default, bool is_global)
  : Statement(pstate), variable_(var), value_(val), is_default_(is_default), is_global_(is_global)
  { statement_type(ASSIGNMENT); }
  Assignment::Assignment(const Assignment* ptr)
  : Statement(ptr),
    variable_(ptr->variable_),
    value_(ptr->value_),
    is_default_(ptr->is_default_),
    is_global_(ptr->is_global_)
  { statement_type(ASSIGNMENT); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Import::Import(ParserState pstate)
  : Statement(pstate),
    urls_(std::vector<Expression_Obj>()),
    incs_(std::vector<Include>()),
    import_queries_()
  { statement_type(IMPORT); }
  Import::Import(const Import* ptr)
  : Statement(ptr),
    urls_(ptr->urls_),
    incs_(ptr->incs_),
    import_queries_(ptr->import_queries_)
  { statement_type(IMPORT); }

  std::vector<Include>& Import::incs() { return incs_; }
  std::vector<Expression_Obj>& Import::urls() { return urls_; }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Import_Stub::Import_Stub(ParserState pstate, Include res)
  : Statement(pstate), resource_(res)
  { statement_type(IMPORT_STUB); }
  Import_Stub::Import_Stub(const Import_Stub* ptr)
  : Statement(ptr), resource_(ptr->resource_)
  { statement_type(IMPORT_STUB); }
  Include Import_Stub::resource() { return resource_; };
  std::string Import_Stub::imp_path() { return resource_.imp_path; };
  std::string Import_Stub::abs_path() { return resource_.abs_path; };

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Warning::Warning(ParserState pstate, Expression_Obj msg)
  : Statement(pstate), message_(msg)
  { statement_type(WARNING); }
  Warning::Warning(const Warning* ptr)
  : Statement(ptr), message_(ptr->message_)
  { statement_type(WARNING); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Error::Error(ParserState pstate, Expression_Obj msg)
  : Statement(pstate), message_(msg)
  { statement_type(ERROR); }
  Error::Error(const Error* ptr)
  : Statement(ptr), message_(ptr->message_)
  { statement_type(ERROR); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Debug::Debug(ParserState pstate, Expression_Obj val)
  : Statement(pstate), value_(val)
  { statement_type(DEBUGSTMT); }
  Debug::Debug(const Debug* ptr)
  : Statement(ptr), value_(ptr->value_)
  { statement_type(DEBUGSTMT); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Comment::Comment(ParserState pstate, String_Obj txt, bool is_important)
  : Statement(pstate), text_(txt), is_important_(is_important)
  { statement_type(COMMENT); }
  Comment::Comment(const Comment* ptr)
  : Statement(ptr),
    text_(ptr->text_),
    is_important_(ptr->is_important_)
  { statement_type(COMMENT); }

  bool Comment::is_invisible() const
  {
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  If::If(ParserState pstate, Expression_Obj pred, Block_Obj con, Block_Obj alt)
  : Has_Block(pstate, con), predicate_(pred), alternative_(alt)
  { statement_type(IF); }
  If::If(const If* ptr)
  : Has_Block(ptr),
    predicate_(ptr->predicate_),
    alternative_(ptr->alternative_)
  { statement_type(IF); }

  bool If::has_content()
  {
    return Has_Block::has_content() || (alternative_ && alternative_->has_content());
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  For::For(ParserState pstate,
      std::string var, Expression_Obj lo, Expression_Obj hi, Block_Obj b, bool inc)
  : Has_Block(pstate, b),
    variable_(var), lower_bound_(lo), upper_bound_(hi), is_inclusive_(inc)
  { statement_type(FOR); }
  For::For(const For* ptr)
  : Has_Block(ptr),
    variable_(ptr->variable_),
    lower_bound_(ptr->lower_bound_),
    upper_bound_(ptr->upper_bound_),
    is_inclusive_(ptr->is_inclusive_)
  { statement_type(FOR); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Each::Each(ParserState pstate, std::vector<std::string> vars, Expression_Obj lst, Block_Obj b)
  : Has_Block(pstate, b), variables_(vars), list_(lst)
  { statement_type(EACH); }
  Each::Each(const Each* ptr)
  : Has_Block(ptr), variables_(ptr->variables_), list_(ptr->list_)
  { statement_type(EACH); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  While::While(ParserState pstate, Expression_Obj pred, Block_Obj b)
  : Has_Block(pstate, b), predicate_(pred)
  { statement_type(WHILE); }
  While::While(const While* ptr)
  : Has_Block(ptr), predicate_(ptr->predicate_)
  { statement_type(WHILE); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Return::Return(ParserState pstate, Expression_Obj val)
  : Statement(pstate), value_(val)
  { statement_type(RETURN); }
  Return::Return(const Return* ptr)
  : Statement(ptr), value_(ptr->value_)
  { statement_type(RETURN); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Extension::Extension(ParserState pstate, Selector_List_Obj s)
  : Statement(pstate), selector_(s)
  { statement_type(EXTEND); }
  Extension::Extension(const Extension* ptr)
  : Statement(ptr), selector_(ptr->selector_)
  { statement_type(EXTEND); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Definition::Definition(const Definition* ptr)
  : Has_Block(ptr),
    name_(ptr->name_),
    parameters_(ptr->parameters_),
    environment_(ptr->environment_),
    type_(ptr->type_),
    native_function_(ptr->native_function_),
    c_function_(ptr->c_function_),
    cookie_(ptr->cookie_),
    is_overload_stub_(ptr->is_overload_stub_),
    signature_(ptr->signature_)
  { }

  Definition::Definition(ParserState pstate,
              std::string n,
              Parameters_Obj params,
              Block_Obj b,
              Type t)
  : Has_Block(pstate, b),
    name_(n),
    parameters_(params),
    environment_(0),
    type_(t),
    native_function_(0),
    c_function_(0),
    cookie_(0),
    is_overload_stub_(false),
    signature_(0)
  { }

  Definition::Definition(ParserState pstate,
              Signature sig,
              std::string n,
              Parameters_Obj params,
              Native_Function func_ptr,
              bool overload_stub)
  : Has_Block(pstate, {}),
    name_(n),
    parameters_(params),
    environment_(0),
    type_(FUNCTION),
    native_function_(func_ptr),
    c_function_(0),
    cookie_(0),
    is_overload_stub_(overload_stub),
    signature_(sig)
  { }

  Definition::Definition(ParserState pstate,
              Signature sig,
              std::string n,
              Parameters_Obj params,
              Sass_Function_Entry c_func)
  : Has_Block(pstate, {}),
    name_(n),
    parameters_(params),
    environment_(0),
    type_(FUNCTION),
    native_function_(0),
    c_function_(c_func),
    cookie_(sass_function_get_cookie(c_func)),
    is_overload_stub_(false),
    signature_(sig)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Mixin_Call::Mixin_Call(ParserState pstate, std::string n, Arguments_Obj args, Parameters_Obj b_params, Block_Obj b)
  : Has_Block(pstate, b), name_(n), arguments_(args), block_parameters_(b_params)
  { }
  Mixin_Call::Mixin_Call(const Mixin_Call* ptr)
  : Has_Block(ptr),
    name_(ptr->name_),
    arguments_(ptr->arguments_),
    block_parameters_(ptr->block_parameters_)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Content::Content(ParserState pstate, Arguments_Obj args)
  : Statement(pstate),
    arguments_(args)
  { statement_type(CONTENT); }
  Content::Content(const Content* ptr)
  : Statement(ptr),
    arguments_(ptr->arguments_)
  { statement_type(CONTENT); }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Expression::Expression(ParserState pstate, bool d, bool e, bool i, Type ct)
  : AST_Node(pstate),
    is_delayed_(d),
    is_expanded_(e),
    is_interpolant_(i),
    concrete_type_(ct)
  { }

  Expression::Expression(const Expression* ptr)
  : AST_Node(ptr),
    is_delayed_(ptr->is_delayed_),
    is_expanded_(ptr->is_expanded_),
    is_interpolant_(ptr->is_interpolant_),
    concrete_type_(ptr->concrete_type_)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Unary_Expression::Unary_Expression(ParserState pstate, Type t, Expression_Obj o)
  : Expression(pstate), optype_(t), operand_(o), hash_(0)
  { }
  Unary_Expression::Unary_Expression(const Unary_Expression* ptr)
  : Expression(ptr),
    optype_(ptr->optype_),
    operand_(ptr->operand_),
    hash_(ptr->hash_)
  { }
  const std::string Unary_Expression::type_name() {
    switch (optype_) {
      case PLUS: return "plus";
      case MINUS: return "minus";
      case SLASH: return "slash";
      case NOT: return "not";
      default: return "invalid";
    }
  }
  bool Unary_Expression::operator==(const Expression& rhs) const
  {
    try
    {
      Unary_Expression_Ptr_Const m = Cast<Unary_Expression>(&rhs);
      if (m == 0) return false;
      return type() == m->type() &&
              *operand() == *m->operand();
    }
    catch (std::bad_cast&)
    {
      return false;
    }
    catch (...) { throw; }
  }
  size_t Unary_Expression::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<size_t>()(optype_);
      hash_combine(hash_, operand()->hash());
    };
    return hash_;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Argument::Argument(ParserState pstate, Expression_Obj val, std::string n, bool rest, bool keyword)
  : Expression(pstate), value_(val), name_(n), is_rest_argument_(rest), is_keyword_argument_(keyword), hash_(0)
  {
    if (!name_.empty() && is_rest_argument_) {
      coreError("variable-length argument may not be passed by name", pstate_);
    }
  }
  Argument::Argument(const Argument* ptr)
  : Expression(ptr),
    value_(ptr->value_),
    name_(ptr->name_),
    is_rest_argument_(ptr->is_rest_argument_),
    is_keyword_argument_(ptr->is_keyword_argument_),
    hash_(ptr->hash_)
  {
    if (!name_.empty() && is_rest_argument_) {
      coreError("variable-length argument may not be passed by name", pstate_);
    }
  }

  void Argument::set_delayed(bool delayed)
  {
    if (value_) value_->set_delayed(delayed);
    is_delayed(delayed);
  }

  bool Argument::operator==(const Expression& rhs) const
  {
    try
    {
      Argument_Ptr_Const m = Cast<Argument>(&rhs);
      if (!(m && name() == m->name())) return false;
      return *value() == *m->value();
    }
    catch (std::bad_cast&)
    {
      return false;
    }
    catch (...) { throw; }
  }

  size_t Argument::hash() const
  {
    if (hash_ == 0) {
      hash_ = std::hash<std::string>()(name());
      hash_combine(hash_, value()->hash());
    }
    return hash_;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Arguments::Arguments(ParserState pstate)
  : Expression(pstate),
    Vectorized<Argument_Obj>(),
    has_named_arguments_(false),
    has_rest_argument_(false),
    has_keyword_argument_(false)
  { }
  Arguments::Arguments(const Arguments* ptr)
  : Expression(ptr),
    Vectorized<Argument_Obj>(*ptr),
    has_named_arguments_(ptr->has_named_arguments_),
    has_rest_argument_(ptr->has_rest_argument_),
    has_keyword_argument_(ptr->has_keyword_argument_)
  { }

  void Arguments::set_delayed(bool delayed)
  {
    for (Argument_Obj arg : elements()) {
      if (arg) arg->set_delayed(delayed);
    }
    is_delayed(delayed);
  }

  Argument_Obj Arguments::get_rest_argument()
  {
    if (this->has_rest_argument()) {
      for (Argument_Obj arg : this->elements()) {
        if (arg->is_rest_argument()) {
          return arg;
        }
      }
    }
    return {};
  }

  Argument_Obj Arguments::get_keyword_argument()
  {
    if (this->has_keyword_argument()) {
      for (Argument_Obj arg : this->elements()) {
        if (arg->is_keyword_argument()) {
          return arg;
        }
      }
    }
    return {};
  }

  void Arguments::adjust_after_pushing(Argument_Obj a)
  {
    if (!a->name().empty()) {
      if (has_keyword_argument()) {
        coreError("named arguments must precede variable-length argument", a->pstate());
      }
      has_named_arguments(true);
    }
    else if (a->is_rest_argument()) {
      if (has_rest_argument()) {
        coreError("functions and mixins may only be called with one variable-length argument", a->pstate());
      }
      if (has_keyword_argument_) {
        coreError("only keyword arguments may follow variable arguments", a->pstate());
      }
      has_rest_argument(true);
    }
    else if (a->is_keyword_argument()) {
      if (has_keyword_argument()) {
        coreError("functions and mixins may only be called with one keyword argument", a->pstate());
      }
      has_keyword_argument(true);
    }
    else {
      if (has_rest_argument()) {
        coreError("ordinal arguments must precede variable-length arguments", a->pstate());
      }
      if (has_named_arguments()) {
        coreError("ordinal arguments must precede named arguments", a->pstate());
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Media_Query::Media_Query(ParserState pstate, String_Obj t, size_t s, bool n, bool r)
  : Expression(pstate), Vectorized<Media_Query_Expression_Obj>(s),
    media_type_(t), is_negated_(n), is_restricted_(r)
  { }
  Media_Query::Media_Query(const Media_Query* ptr)
  : Expression(ptr),
    Vectorized<Media_Query_Expression_Obj>(*ptr),
    media_type_(ptr->media_type_),
    is_negated_(ptr->is_negated_),
    is_restricted_(ptr->is_restricted_)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Media_Query_Expression::Media_Query_Expression(ParserState pstate,
                          Expression_Obj f, Expression_Obj v, bool i)
  : Expression(pstate), feature_(f), value_(v), is_interpolated_(i)
  { }
  Media_Query_Expression::Media_Query_Expression(const Media_Query_Expression* ptr)
  : Expression(ptr),
    feature_(ptr->feature_),
    value_(ptr->value_),
    is_interpolated_(ptr->is_interpolated_)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  At_Root_Query::At_Root_Query(ParserState pstate, Expression_Obj f, Expression_Obj v, bool i)
  : Expression(pstate), feature_(f), value_(v)
  { }
  At_Root_Query::At_Root_Query(const At_Root_Query* ptr)
  : Expression(ptr),
    feature_(ptr->feature_),
    value_(ptr->value_)
  { }

  bool At_Root_Query::exclude(std::string str)
  {
    bool with = feature() && unquote(feature()->to_string()).compare("with") == 0;
    List_Ptr l = static_cast<List_Ptr>(value().ptr());
    std::string v;

    if (with)
    {
      if (!l || l->length() == 0) return str.compare("rule") != 0;
      for (size_t i = 0, L = l->length(); i < L; ++i)
      {
        v = unquote((*l)[i]->to_string());
        if (v.compare("all") == 0 || v == str) return false;
      }
      return true;
    }
    else
    {
      if (!l || !l->length()) return str.compare("rule") == 0;
      for (size_t i = 0, L = l->length(); i < L; ++i)
      {
        v = unquote((*l)[i]->to_string());
        if (v.compare("all") == 0 || v == str) return true;
      }
      return false;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  At_Root_Block::At_Root_Block(ParserState pstate, Block_Obj b, At_Root_Query_Obj e)
  : Has_Block(pstate, b), expression_(e)
  { statement_type(ATROOT); }
  At_Root_Block::At_Root_Block(const At_Root_Block* ptr)
  : Has_Block(ptr), expression_(ptr->expression_)
  { statement_type(ATROOT); }

  bool At_Root_Block::bubbles() {
    return true;
  }

  bool At_Root_Block::exclude_node(Statement_Obj s) {
    if (expression() == 0)
    {
      return s->statement_type() == Statement::RULESET;
    }

    if (s->statement_type() == Statement::DIRECTIVE)
    {
      if (Directive_Obj dir = Cast<Directive>(s))
      {
        std::string keyword(dir->keyword());
        if (keyword.length() > 0) keyword.erase(0, 1);
        return expression()->exclude(keyword);
      }
    }
    if (s->statement_type() == Statement::MEDIA)
    {
      return expression()->exclude("media");
    }
    if (s->statement_type() == Statement::RULESET)
    {
      return expression()->exclude("rule");
    }
    if (s->statement_type() == Statement::SUPPORTS)
    {
      return expression()->exclude("supports");
    }
    if (Directive_Obj dir = Cast<Directive>(s))
    {
      if (dir->is_keyframes()) return expression()->exclude("keyframes");
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Parameter::Parameter(ParserState pstate, std::string n, Expression_Obj def, bool rest)
  : AST_Node(pstate), name_(n), default_value_(def), is_rest_parameter_(rest)
  { }
  Parameter::Parameter(const Parameter* ptr)
  : AST_Node(ptr),
    name_(ptr->name_),
    default_value_(ptr->default_value_),
    is_rest_parameter_(ptr->is_rest_parameter_)
  { }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Parameters::Parameters(ParserState pstate)
  : AST_Node(pstate),
    Vectorized<Parameter_Obj>(),
    has_optional_parameters_(false),
    has_rest_parameter_(false)
  { }
  Parameters::Parameters(const Parameters* ptr)
  : AST_Node(ptr),
    Vectorized<Parameter_Obj>(*ptr),
    has_optional_parameters_(ptr->has_optional_parameters_),
    has_rest_parameter_(ptr->has_rest_parameter_)
  { }

  void Parameters::adjust_after_pushing(Parameter_Obj p)
  {
    if (p->default_value()) {
      if (has_rest_parameter()) {
        coreError("optional parameters may not be combined with variable-length parameters", p->pstate());
      }
      has_optional_parameters(true);
    }
    else if (p->is_rest_parameter()) {
      if (has_rest_parameter()) {
        coreError("functions and mixins cannot have more than one variable-length parameter", p->pstate());
      }
      has_rest_parameter(true);
    }
    else {
      if (has_rest_parameter()) {
        coreError("required parameters must precede variable-length parameters", p->pstate());
      }
      if (has_optional_parameters()) {
        coreError("required parameters must precede optional parameters", p->pstate());
      }
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  IMPLEMENT_AST_OPERATORS(Ruleset);
  IMPLEMENT_AST_OPERATORS(Media_Block);
  IMPLEMENT_AST_OPERATORS(Import);
  IMPLEMENT_AST_OPERATORS(Import_Stub);
  IMPLEMENT_AST_OPERATORS(Directive);
  IMPLEMENT_AST_OPERATORS(At_Root_Block);
  IMPLEMENT_AST_OPERATORS(While);
  IMPLEMENT_AST_OPERATORS(Each);
  IMPLEMENT_AST_OPERATORS(For);
  IMPLEMENT_AST_OPERATORS(If);
  IMPLEMENT_AST_OPERATORS(Mixin_Call);
  IMPLEMENT_AST_OPERATORS(Extension);
  IMPLEMENT_AST_OPERATORS(Media_Query);
  IMPLEMENT_AST_OPERATORS(Media_Query_Expression);
  IMPLEMENT_AST_OPERATORS(Debug);
  IMPLEMENT_AST_OPERATORS(Error);
  IMPLEMENT_AST_OPERATORS(Warning);
  IMPLEMENT_AST_OPERATORS(Assignment);
  IMPLEMENT_AST_OPERATORS(Return);
  IMPLEMENT_AST_OPERATORS(At_Root_Query);
  IMPLEMENT_AST_OPERATORS(Comment);
  IMPLEMENT_AST_OPERATORS(Parameters);
  IMPLEMENT_AST_OPERATORS(Parameter);
  IMPLEMENT_AST_OPERATORS(Arguments);
  IMPLEMENT_AST_OPERATORS(Argument);
  IMPLEMENT_AST_OPERATORS(Unary_Expression);
  IMPLEMENT_AST_OPERATORS(Block);
  IMPLEMENT_AST_OPERATORS(Content);
  IMPLEMENT_AST_OPERATORS(Trace);
  IMPLEMENT_AST_OPERATORS(Keyframe_Rule);
  IMPLEMENT_AST_OPERATORS(Bubble);
  IMPLEMENT_AST_OPERATORS(Definition);
  IMPLEMENT_AST_OPERATORS(Declaration);

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  CheckNesting::CheckNesting()
  : parents(std::vector<Statement_Ptr>()),
    traces(std::vector<Backtrace>()),
    parent(0), current_mixin_definition(0)
  { }

  void error(AST_Node_Ptr node, Backtraces traces, std::string msg) {
    traces.push_back(Backtrace(node->pstate()));
    throw Exception::InvalidSass(node->pstate(), traces, msg);
  }

  Statement_Ptr CheckNesting::visit_children(Statement_Ptr parent)
  {
    Statement_Ptr old_parent = this->parent;

    if (At_Root_Block_Ptr root = Cast<At_Root_Block>(parent)) {
      std::vector<Statement_Ptr> old_parents = this->parents;
      std::vector<Statement_Ptr> new_parents;

      for (size_t i = 0, L = this->parents.size(); i < L; i++) {
        Statement_Ptr p = this->parents.at(i);
        if (!root->exclude_node(p)) {
          new_parents.push_back(p);
        }
      }
      this->parents = new_parents;

      for (size_t i = this->parents.size(); i > 0; i--) {
        Statement_Ptr p = 0;
        Statement_Ptr gp = 0;
        if (i > 0) p = this->parents.at(i - 1);
        if (i > 1) gp = this->parents.at(i - 2);

        if (!this->is_transparent_parent(p, gp)) {
          this->parent = p;
          break;
        }
      }

      At_Root_Block_Ptr ar = Cast<At_Root_Block>(parent);
      Block_Ptr ret = ar->block();

      if (ret != NULL) {
        for (auto n : ret->elements()) {
          n->perform(this);
        }
      }

      this->parent = old_parent;
      this->parents = old_parents;

      return ret;
    }

    if (!this->is_transparent_parent(parent, old_parent)) {
      this->parent = parent;
    }

    this->parents.push_back(parent);

    Block_Ptr b = Cast<Block>(parent);

    if (Trace_Ptr trace = Cast<Trace>(parent)) {
      if (trace->type() == 'i') {
        this->traces.push_back(Backtrace(trace->pstate()));
      }
    }

    if (!b) {
      if (Has_Block_Ptr bb = Cast<Has_Block>(parent)) {
        b = bb->block();
      }
    }

    if (b) {
      for (auto n : b->elements()) {
        n->perform(this);
      }
    }

    this->parent = old_parent;
    this->parents.pop_back();

    if (Trace_Ptr trace = Cast<Trace>(parent)) {
      if (trace->type() == 'i') {
        this->traces.pop_back();
      }
    }

    return b;
  }


  Statement_Ptr CheckNesting::operator()(Block_Ptr b)
  {
    return this->visit_children(b);
  }

  Statement_Ptr CheckNesting::operator()(Definition_Ptr n)
  {
    if (!this->should_visit(n)) return NULL;
    if (!is_mixin(n)) {
      visit_children(n);
      return n;
    }

    Definition_Ptr old_mixin_definition = this->current_mixin_definition;
    this->current_mixin_definition = n;

    visit_children(n);

    this->current_mixin_definition = old_mixin_definition;

    return n;
  }

  Statement_Ptr CheckNesting::operator()(If_Ptr i)
  {
    this->visit_children(i);

    if (Block_Ptr b = Cast<Block>(i->alternative())) {
      for (auto n : b->elements()) n->perform(this);
    }

    return i;
  }

  bool CheckNesting::should_visit(Statement_Ptr node)
  {
    if (!this->parent) return true;

    if (Cast<Content>(node))
    { this->invalid_content_parent(this->parent, node); }

    if (is_charset(node))
    { this->invalid_charset_parent(this->parent, node); }

    if (Cast<Extension>(node))
    { this->invalid_extend_parent(this->parent, node); }

    // if (Cast<Import>(node))
    // { this->invalid_import_parent(this->parent); }

    if (this->is_mixin(node))
    { this->invalid_mixin_definition_parent(this->parent, node); }

    if (this->is_function(node))
    { this->invalid_function_parent(this->parent, node); }

    if (this->is_function(this->parent))
    { this->invalid_function_child(node); }

    if (Declaration_Ptr d = Cast<Declaration>(node))
    {
      this->invalid_prop_parent(this->parent, node);
      this->invalid_value_child(d->value());
    }

    if (Cast<Declaration>(this->parent))
    { this->invalid_prop_child(node); }

    if (Cast<Return>(node))
    { this->invalid_return_parent(this->parent, node); }

    return true;
  }

  void CheckNesting::invalid_content_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    if (!this->current_mixin_definition) {
      error(node, traces, "@content may only be used within a mixin.");
    }
  }

  void CheckNesting::invalid_charset_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    if (!(
        is_root_node(parent)
    )) {
      error(node, traces, "@charset may only be used at the root of a document.");
    }
  }

  void CheckNesting::invalid_extend_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    if (!(
        Cast<Ruleset>(parent) ||
        Cast<Mixin_Call>(parent) ||
        is_mixin(parent)
    )) {
      error(node, traces, "Extend directives may only be used within rules.");
    }
  }

  // void CheckNesting::invalid_import_parent(Statement_Ptr parent, AST_Node_Ptr node)
  // {
  //   for (auto pp : this->parents) {
  //     if (
  //         Cast<Each>(pp) ||
  //         Cast<For>(pp) ||
  //         Cast<If>(pp) ||
  //         Cast<While>(pp) ||
  //         Cast<Trace>(pp) ||
  //         Cast<Mixin_Call>(pp) ||
  //         is_mixin(pp)
  //     ) {
  //       error(node, traces, "Import directives may not be defined within control directives or other mixins.");
  //     }
  //   }

  //   if (this->is_root_node(parent)) {
  //     return;
  //   }

  //   if (false/*n.css_import?*/) {
  //     error(node, traces, "CSS import directives may only be used at the root of a document.");
  //   }
  // }

  void CheckNesting::invalid_mixin_definition_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    for (Statement_Ptr pp : this->parents) {
      if (
          Cast<Each>(pp) ||
          Cast<For>(pp) ||
          Cast<If>(pp) ||
          Cast<While>(pp) ||
          Cast<Trace>(pp) ||
          Cast<Mixin_Call>(pp) ||
          is_mixin(pp)
      ) {
        error(node, traces, "Mixins may not be defined within control directives or other mixins.");
      }
    }
  }

  void CheckNesting::invalid_function_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    for (Statement_Ptr pp : this->parents) {
      if (
          Cast<Each>(pp) ||
          Cast<For>(pp) ||
          Cast<If>(pp) ||
          Cast<While>(pp) ||
          Cast<Trace>(pp) ||
          Cast<Mixin_Call>(pp) ||
          is_mixin(pp)
      ) {
        error(node, traces, "Functions may not be defined within control directives or other mixins.");
      }
    }
  }

  void CheckNesting::invalid_function_child(Statement_Ptr child)
  {
    if (!(
        Cast<Each>(child) ||
        Cast<For>(child) ||
        Cast<If>(child) ||
        Cast<While>(child) ||
        Cast<Trace>(child) ||
        Cast<Comment>(child) ||
        Cast<Debug>(child) ||
        Cast<Return>(child) ||
        Cast<Variable>(child) ||
        // Ruby Sass doesn't distinguish variables and assignments
        Cast<Assignment>(child) ||
        Cast<Warning>(child) ||
        Cast<Error>(child)
    )) {
      error(child, traces, "Functions can only contain variable declarations and control directives.");
    }
  }

  void CheckNesting::invalid_prop_child(Statement_Ptr child)
  {
    if (!(
        Cast<Each>(child) ||
        Cast<For>(child) ||
        Cast<If>(child) ||
        Cast<While>(child) ||
        Cast<Trace>(child) ||
        Cast<Comment>(child) ||
        Cast<Declaration>(child) ||
        Cast<Mixin_Call>(child)
    )) {
      error(child, traces, "Illegal nesting: Only properties may be nested beneath properties.");
    }
  }

  void CheckNesting::invalid_prop_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    if (!(
        is_mixin(parent) ||
        is_directive_node(parent) ||
        Cast<Ruleset>(parent) ||
        Cast<Keyframe_Rule>(parent) ||
        Cast<Declaration>(parent) ||
        Cast<Mixin_Call>(parent)
    )) {
      error(node, traces, "Properties are only allowed within rules, directives, mixin includes, or other properties.");
    }
  }

  void CheckNesting::invalid_value_child(AST_Node_Ptr d)
  {
    if (Map_Ptr m = Cast<Map>(d)) {
      traces.push_back(Backtrace(m->pstate()));
      throw Exception::InvalidValue(traces, *m);
    }
    if (Number_Ptr n = Cast<Number>(d)) {
      if (!n->is_valid_css_unit()) {
        traces.push_back(Backtrace(n->pstate()));
        throw Exception::InvalidValue(traces, *n);
      }
    }

    // error(dbg + " isn't a valid CSS value.", m->pstate(),);

  }

  void CheckNesting::invalid_return_parent(Statement_Ptr parent, AST_Node_Ptr node)
  {
    if (!this->is_function(parent)) {
      error(node, traces, "@return may only be used within a function.");
    }
  }

  bool CheckNesting::is_transparent_parent(Statement_Ptr parent, Statement_Ptr grandparent)
  {
    bool parent_bubbles = parent && parent->bubbles();

    bool valid_bubble_node = parent_bubbles &&
                             !is_root_node(grandparent) &&
                             !is_at_root_node(grandparent);

    return Cast<Import>(parent) ||
           Cast<Each>(parent) ||
           Cast<For>(parent) ||
           Cast<If>(parent) ||
           Cast<While>(parent) ||
           Cast<Trace>(parent) ||
           valid_bubble_node;
  }

  bool CheckNesting::is_charset(Statement_Ptr n)
  {
    Directive_Ptr d = Cast<Directive>(n);
    return d && d->keyword() == "charset";
  }

  bool CheckNesting::is_mixin(Statement_Ptr n)
  {
    Definition_Ptr def = Cast<Definition>(n);
    return def && def->type() == Definition::MIXIN;
  }

  bool CheckNesting::is_function(Statement_Ptr n)
  {
    Definition_Ptr def = Cast<Definition>(n);
    return def && def->type() == Definition::FUNCTION;
  }

  bool CheckNesting::is_root_node(Statement_Ptr n)
  {
    if (Cast<Ruleset>(n)) return false;

    Block_Ptr b = Cast<Block>(n);
    return b && b->is_root();
  }

  bool CheckNesting::is_at_root_node(Statement_Ptr n)
  {
    return Cast<At_Root_Block>(n) != NULL;
  }

  bool CheckNesting::is_directive_node(Statement_Ptr n)
  {
    return Cast<Directive>(n) ||
           Cast<Import>(n) ||
           Cast<Media_Block>(n) ||
           Cast<Supports_Block>(n);
  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {


  Node Node::createCombinator(const Complex_Selector::Combinator& combinator) {
    NodeDequePtr null;
    return Node(COMBINATOR, combinator, NULL /*pSelector*/, null /*pCollection*/);
  }


  Node Node::createSelector(const Complex_Selector& pSelector) {
    NodeDequePtr null;

    Complex_Selector_Ptr pStripped = SASS_MEMORY_COPY(&pSelector);
    pStripped->tail({});
    pStripped->combinator(Complex_Selector::ANCESTOR_OF);

    Node n(SELECTOR, Complex_Selector::ANCESTOR_OF, pStripped, null /*pCollection*/);
    n.got_line_feed = pSelector.has_line_feed();
    return n;
  }


  Node Node::createCollection() {
    NodeDequePtr pEmptyCollection = std::make_shared<NodeDeque>();
    return Node(COLLECTION, Complex_Selector::ANCESTOR_OF, NULL /*pSelector*/, pEmptyCollection);
  }


  Node Node::createCollection(const NodeDeque& values) {
    NodeDequePtr pShallowCopiedCollection = std::make_shared<NodeDeque>(values);
    return Node(COLLECTION, Complex_Selector::ANCESTOR_OF, NULL /*pSelector*/, pShallowCopiedCollection);
  }


  Node Node::createNil() {
    NodeDequePtr null;
    return Node(NIL, Complex_Selector::ANCESTOR_OF, NULL /*pSelector*/, null /*pCollection*/);
  }


  Node::Node(const TYPE& type, Complex_Selector::Combinator combinator, Complex_Selector_Ptr pSelector, NodeDequePtr& pCollection)
  : got_line_feed(false), mType(type), mCombinator(combinator), mpSelector(pSelector), mpCollection(pCollection)
  { if (pSelector) got_line_feed = pSelector->has_line_feed(); }


  Node Node::klone() const {
    NodeDequePtr pNewCollection = std::make_shared<NodeDeque>();
    if (mpCollection) {
      for (NodeDeque::iterator iter = mpCollection->begin(), iterEnd = mpCollection->end(); iter != iterEnd; iter++) {
        Node& toClone = *iter;
        pNewCollection->push_back(toClone.klone());
      }
    }

    Node n(mType, mCombinator, mpSelector ? SASS_MEMORY_COPY(mpSelector) : NULL, pNewCollection);
    n.got_line_feed = got_line_feed;
    return n;
  }


  bool Node::contains(const Node& potentialChild) const {
    bool found = false;

    for (NodeDeque::iterator iter = mpCollection->begin(), iterEnd = mpCollection->end(); iter != iterEnd; iter++) {
      Node& toTest = *iter;

      if (toTest == potentialChild) {
        found = true;
        break;
      }
    }

    return found;
  }


  bool Node::operator==(const Node& rhs) const {
    if (this->type() != rhs.type()) {
      return false;
    }

    if (this->isCombinator()) {

      return this->combinator() == rhs.combinator();

    } else if (this->isNil()) {

      return true; // no state to check

    } else if (this->isSelector()){

      return *this->selector() == *rhs.selector();

    } else if (this->isCollection()) {

      if (this->collection()->size() != rhs.collection()->size()) {
        return false;
      }

      for (NodeDeque::iterator lhsIter = this->collection()->begin(), lhsIterEnd = this->collection()->end(),
           rhsIter = rhs.collection()->begin(); lhsIter != lhsIterEnd; lhsIter++, rhsIter++) {

        if (*lhsIter != *rhsIter) {
          return false;
        }

      }

      return true;

    }

    // We shouldn't get here.
    throw "Comparing unknown node types. A new type was probably added and this method wasn't implemented for it.";
  }


  void Node::plus(Node& rhs) {
    if (!this->isCollection() || !rhs.isCollection()) {
      throw "Both the current node and rhs must be collections.";
    }
    this->collection()->insert(this->collection()->end(), rhs.collection()->begin(), rhs.collection()->end());
  }

#ifdef DEBUG
  std::ostream& operator<<(std::ostream& os, const Node& node) {

    if (node.isCombinator()) {

      switch (node.combinator()) {
        case Complex_Selector::ANCESTOR_OF: os << "\" \""; break;
        case Complex_Selector::PARENT_OF:   os << "\">\""; break;
        case Complex_Selector::PRECEDES:    os << "\"~\""; break;
        case Complex_Selector::ADJACENT_TO: os << "\"+\""; break;
        case Complex_Selector::REFERENCE: os    << "\"/\""; break;
      }

    } else if (node.isNil()) {

      os << "nil";

    } else if (node.isSelector()){

      os << node.selector()->head()->to_string();

    } else if (node.isCollection()) {

      os << "[";

      for (NodeDeque::iterator iter = node.collection()->begin(), iterBegin = node.collection()->begin(), iterEnd = node.collection()->end(); iter != iterEnd; iter++) {
        if (iter != iterBegin) {
          os << ", ";
        }

        os << (*iter);
      }

      os << "]";

    }

    return os;

  }
#endif


  Node complexSelectorToNode(Complex_Selector_Ptr pToConvert) {
    if (pToConvert == NULL) {
      return Node::createNil();
    }
    Node node = Node::createCollection();
    node.got_line_feed = pToConvert->has_line_feed();
    bool has_lf = pToConvert->has_line_feed();

    // unwrap the selector from parent ref
    if (pToConvert->head() && pToConvert->head()->has_parent_ref()) {
      Complex_Selector_Obj tail = pToConvert->tail();
      if (tail) tail->has_line_feed(pToConvert->has_line_feed());
      pToConvert = tail;
    }

    while (pToConvert) {

      bool empty_parent_ref = pToConvert->head() && pToConvert->head()->is_empty_reference();

      // the first Complex_Selector may contain a dummy head pointer, skip it.
      if (pToConvert->head() && !empty_parent_ref) {
        node.collection()->push_back(Node::createSelector(*pToConvert));
        if (has_lf) node.collection()->back().got_line_feed = has_lf;
        if (pToConvert->head() || empty_parent_ref) {
          if (pToConvert->tail()) {
            pToConvert->tail()->has_line_feed(pToConvert->has_line_feed());
          }
        }
        has_lf = false;
      }

      if (pToConvert->combinator() != Complex_Selector::ANCESTOR_OF) {
        node.collection()->push_back(Node::createCombinator(pToConvert->combinator()));
        if (has_lf) node.collection()->back().got_line_feed = has_lf;
        has_lf = false;
      }

      if (pToConvert && empty_parent_ref && pToConvert->tail()) {
        // pToConvert->tail()->has_line_feed(pToConvert->has_line_feed());
      }

      pToConvert = pToConvert->tail();
    }

    return node;
  }


  Complex_Selector_Ptr nodeToComplexSelector(const Node& toConvert) {
    if (toConvert.isNil()) {
      return NULL;
    }


    if (!toConvert.isCollection()) {
      throw "The node to convert to a Complex_Selector_Ptr must be a collection type or nil.";
    }


    NodeDeque& childNodes = *toConvert.collection();

    std::string noPath("");
    Complex_Selector_Obj pFirst = SASS_MEMORY_NEW(Complex_Selector, ParserState("[NODE]"), Complex_Selector::ANCESTOR_OF, {}, {});

    Complex_Selector_Obj pCurrent = pFirst;

    if (toConvert.isSelector()) pFirst->has_line_feed(toConvert.got_line_feed);
    if (toConvert.isCombinator()) pFirst->has_line_feed(toConvert.got_line_feed);

    for (NodeDeque::iterator childIter = childNodes.begin(), childIterEnd = childNodes.end(); childIter != childIterEnd; childIter++) {

      Node& child = *childIter;

      if (child.isSelector()) {
        // JMA - need to clone the selector, because they can end up getting shared across Node
        // collections, and can result in an infinite loop during the call to parentSuperselector()
        pCurrent->tail(SASS_MEMORY_COPY(child.selector()));
        // if (child.got_line_feed) pCurrent->has_line_feed(child.got_line_feed);
        pCurrent = pCurrent->tail();
      } else if (child.isCombinator()) {
        pCurrent->combinator(child.combinator());
        if (child.got_line_feed) pCurrent->has_line_feed(child.got_line_feed);

        // if the next node is also a combinator, create another Complex_Selector to hold it so it doesn't replace the current combinator
        if (childIter+1 != childIterEnd) {
          Node& nextNode = *(childIter+1);
          if (nextNode.isCombinator()) {
            pCurrent->tail(SASS_MEMORY_NEW(Complex_Selector, ParserState("[NODE]"), Complex_Selector::ANCESTOR_OF, {}, {}));
            if (nextNode.got_line_feed) pCurrent->tail()->has_line_feed(nextNode.got_line_feed);
            pCurrent = pCurrent->tail();
          }
        }
      } else {
        throw "The node to convert's children must be only combinators or selectors.";
      }
    }

    // Put the dummy Compound_Selector in the first position, for consistency with the rest of libsass
    Compound_Selector_Ptr fakeHead = SASS_MEMORY_NEW(Compound_Selector, ParserState("[NODE]"), 1);
    Parent_Selector_Ptr selectorRef = SASS_MEMORY_NEW(Parent_Selector, ParserState("[NODE]"));
    fakeHead->elements().push_back(selectorRef);
    if (toConvert.got_line_feed) pFirst->has_line_feed(toConvert.got_line_feed);
    // pFirst->has_line_feed(pFirst->has_line_feed() || pFirst->tail()->has_line_feed() || toConvert.got_line_feed);
    pFirst->head(fakeHead);
    return SASS_MEMORY_COPY(pFirst);
  }

  // A very naive trim function, which removes duplicates in a node
  // This is only used in Complex_Selector::unify_with for now, may need modifications to fit other needs
  Node Node::naiveTrim(Node& seqses) {

    std::vector<Node*> res;
    std::vector<Complex_Selector_Obj> known;

    NodeDeque::reverse_iterator seqsesIter = seqses.collection()->rbegin(),
                                seqsesIterEnd = seqses.collection()->rend();

    for (; seqsesIter != seqsesIterEnd; ++seqsesIter)
    {
      Node& seqs1 = *seqsesIter;
      if( seqs1.isSelector() ) {
        Complex_Selector_Obj sel = seqs1.selector();
        std::vector<Complex_Selector_Obj>::iterator it;
        bool found = false;
        for (it = known.begin(); it != known.end(); ++it) {
          if (**it == *sel) { found = true; break; }
        }
        if( !found ) {
          known.push_back(seqs1.selector());
          res.push_back(&seqs1);
        }
      } else {
        res.push_back(&seqs1);
      }
    }

    Node result = Node::createCollection();

    for (size_t i = res.size() - 1; i != std::string::npos; --i) {
      result.collection()->push_back(*res[i]);
    }

    return result;
  }
}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



#define LFEED "\n"

// C++ helper
namespace Sass {
  // see sass_copy_c_string(std::string str)
  static inline JsonNode* json_mkstream(const std::stringstream& stream)
  {
    // hold on to string on stack!
    std::string str(stream.str());
    return json_mkstring(str.c_str());
  }

  static int handle_error(Sass_Context* c_ctx) {
    try {
      throw;
    }
    catch (Exception::Base& e) {
      std::stringstream msg_stream;
      std::string cwd(Sass::File::get_cwd());
      std::string msg_prefix(e.errtype());
      bool got_newline = false;
      msg_stream << msg_prefix << ": ";
      const char* msg = e.what();
      while (msg && *msg) {
        if (*msg == '\r') {
          got_newline = true;
        }
        else if (*msg == '\n') {
          got_newline = true;
        }
        else if (got_newline) {
          msg_stream << std::string(msg_prefix.size() + 2, ' ');
          got_newline = false;
        }
        msg_stream << *msg;
        ++msg;
      }
      if (!got_newline) msg_stream << "\n";

      if (e.traces.empty()) {
        // we normally should have some traces, still here as a fallback
        std::string rel_path(Sass::File::abs2rel(e.pstate.path, cwd, cwd));
        msg_stream << std::string(msg_prefix.size() + 2, ' ');
        msg_stream << " on line " << e.pstate.line + 1 << " of " << rel_path << "\n";
      }
      else {
        std::string rel_path(Sass::File::abs2rel(e.pstate.path, cwd, cwd));
        msg_stream << traces_to_string(e.traces, "        ");
      }

      // now create the code trace (ToDo: maybe have util functions?)
      if (e.pstate.line != std::string::npos &&
          e.pstate.column != std::string::npos &&
          e.pstate.src != nullptr) {
        size_t lines = e.pstate.line;
        // scan through src until target line
        // move line_beg pointer to line start
        const char* line_beg;
        for (line_beg = e.pstate.src; *line_beg != '\0'; ++line_beg) {
          if (lines == 0) break;
          if (*line_beg == '\n') --lines;
        }
        // move line_end before next newline character
        const char* line_end;
        for (line_end = line_beg; *line_end != '\0'; ++line_end) {
          if (*line_end == '\n' || *line_end == '\r') break;
        }
        if (*line_end != '\0') ++line_end;
        size_t line_len = line_end - line_beg;
        size_t move_in = 0; size_t shorten = 0;
        size_t left_chars = 42; size_t max_chars = 76;
        // reported excerpt should not exceed `max_chars` chars
        if (e.pstate.column > line_len) left_chars = e.pstate.column;
        if (e.pstate.column > left_chars) move_in = e.pstate.column - left_chars;
        if (line_len > max_chars + move_in) shorten = line_len - move_in - max_chars;
        utf8::advance(line_beg, move_in, line_end);
        utf8::retreat(line_end, shorten, line_beg);
        std::string sanitized; std::string marker(e.pstate.column - move_in, '-');
        utf8::replace_invalid(line_beg, line_end, std::back_inserter(sanitized));
        msg_stream << ">> " << sanitized << "\n";
        msg_stream << "   " << marker << "^\n";
      }

      JsonNode* json_err = json_mkobject();
      json_append_member(json_err, "status", json_mknumber(1));
      json_append_member(json_err, "file", json_mkstring(e.pstate.path));
      json_append_member(json_err, "line", json_mknumber((double)(e.pstate.line + 1)));
      json_append_member(json_err, "column", json_mknumber((double)(e.pstate.column + 1)));
      json_append_member(json_err, "message", json_mkstring(e.what()));
      json_append_member(json_err, "formatted", json_mkstream(msg_stream));
      try { c_ctx->error_json = json_stringify(json_err, "  "); }
      catch (...) {}
      c_ctx->error_message = sass_copy_string(msg_stream.str());
      c_ctx->error_text = sass_copy_c_string(e.what());
      c_ctx->error_status = 1;
      c_ctx->error_file = sass_copy_c_string(e.pstate.path);
      c_ctx->error_line = e.pstate.line + 1;
      c_ctx->error_column = e.pstate.column + 1;
      c_ctx->error_src = e.pstate.src;
      c_ctx->output_string = 0;
      c_ctx->source_map_string = 0;
      json_delete(json_err);
    }
    catch (std::bad_alloc& ba) {
      std::stringstream msg_stream;
      JsonNode* json_err = json_mkobject();
      msg_stream << "Unable to allocate memory: " << ba.what() << std::endl;
      json_append_member(json_err, "status", json_mknumber(2));
      json_append_member(json_err, "message", json_mkstring(ba.what()));
      json_append_member(json_err, "formatted", json_mkstream(msg_stream));
      try { c_ctx->error_json = json_stringify(json_err, "  "); }
      catch (...) {}
      c_ctx->error_message = sass_copy_string(msg_stream.str());
      c_ctx->error_text = sass_copy_c_string(ba.what());
      c_ctx->error_status = 2;
      c_ctx->output_string = 0;
      c_ctx->source_map_string = 0;
      json_delete(json_err);
    }
    catch (std::exception& e) {
      std::stringstream msg_stream;
      JsonNode* json_err = json_mkobject();
      msg_stream << "Internal Error: " << e.what() << std::endl;
      json_append_member(json_err, "status", json_mknumber(3));
      json_append_member(json_err, "message", json_mkstring(e.what()));
      json_append_member(json_err, "formatted", json_mkstream(msg_stream));
      try { c_ctx->error_json = json_stringify(json_err, "  "); }
      catch (...) {}
      c_ctx->error_message = sass_copy_string(msg_stream.str());
      c_ctx->error_text = sass_copy_c_string(e.what());
      c_ctx->error_status = 3;
      c_ctx->output_string = 0;
      c_ctx->source_map_string = 0;
      json_delete(json_err);
    }
    catch (std::string& e) {
      std::stringstream msg_stream;
      JsonNode* json_err = json_mkobject();
      msg_stream << "Internal Error: " << e << std::endl;
      json_append_member(json_err, "status", json_mknumber(4));
      json_append_member(json_err, "message", json_mkstring(e.c_str()));
      json_append_member(json_err, "formatted", json_mkstream(msg_stream));
      try { c_ctx->error_json = json_stringify(json_err, "  "); }
      catch (...) {}
      c_ctx->error_message = sass_copy_string(msg_stream.str());
      c_ctx->error_text = sass_copy_c_string(e.c_str());
      c_ctx->error_status = 4;
      c_ctx->output_string = 0;
      c_ctx->source_map_string = 0;
      json_delete(json_err);
    }
    catch (const char* e) {
      std::stringstream msg_stream;
      JsonNode* json_err = json_mkobject();
      msg_stream << "Internal Error: " << e << std::endl;
      json_append_member(json_err, "status", json_mknumber(4));
      json_append_member(json_err, "message", json_mkstring(e));
      json_append_member(json_err, "formatted", json_mkstream(msg_stream));
      try { c_ctx->error_json = json_stringify(json_err, "  "); }
      catch (...) {}
      c_ctx->error_message = sass_copy_string(msg_stream.str());
      c_ctx->error_text = sass_copy_c_string(e);
      c_ctx->error_status = 4;
      c_ctx->output_string = 0;
      c_ctx->source_map_string = 0;
      json_delete(json_err);
    }
    catch (...) {
      std::stringstream msg_stream;
      JsonNode* json_err = json_mkobject();
      msg_stream << "Unknown error occurred" << std::endl;
      json_append_member(json_err, "status", json_mknumber(5));
      json_append_member(json_err, "message", json_mkstring("unknown"));
      try { c_ctx->error_json = json_stringify(json_err, "  "); }
      catch (...) {}
      c_ctx->error_message = sass_copy_string(msg_stream.str());
      c_ctx->error_text = sass_copy_c_string("unknown");
      c_ctx->error_status = 5;
      c_ctx->output_string = 0;
      c_ctx->source_map_string = 0;
      json_delete(json_err);
    }
    return c_ctx->error_status;
  }

  // allow one error handler to throw another error
  // this can happen with invalid utf8 and json lib
  static int handle_errors(Sass_Context* c_ctx) {
    try { return handle_error(c_ctx); }
    catch (...) { return handle_error(c_ctx); }
  }

  static Block_Obj sass_parse_block(Sass_Compiler* compiler) throw()
  {

    // assert valid pointer
    if (compiler == 0) return {};
    // The cpp context must be set by now
    Context* cpp_ctx = compiler->cpp_ctx;
    Sass_Context* c_ctx = compiler->c_ctx;
    // We will take care to wire up the rest
    compiler->cpp_ctx->c_compiler = compiler;
    compiler->state = SASS_COMPILER_PARSED;

    try {

      // get input/output path from options
      std::string input_path = safe_str(c_ctx->input_path);
      std::string output_path = safe_str(c_ctx->output_path);

      // maybe skip some entries of included files
      // we do not include stdin for data contexts
      bool skip = c_ctx->type == SASS_CONTEXT_DATA;

      // dispatch parse call
      Block_Obj root(cpp_ctx->parse());
      // abort on errors
      if (!root) return {};

      // skip all prefixed files? (ToDo: check srcmap)
      // IMO source-maps should point to headers already
      // therefore don't skip it for now. re-enable or
      // remove completely once this is tested
      size_t headers = cpp_ctx->head_imports;

      // copy the included files on to the context (dont forget to free later)
      if (copy_strings(cpp_ctx->get_included_files(skip, headers), &c_ctx->included_files) == NULL)
        throw(std::bad_alloc());

      // return parsed block
      return root;

    }
    // pass errors to generic error handler
    catch (...) { handle_errors(c_ctx); }

    // error
    return {};

  }

}

extern "C" {
  using namespace Sass;

  static void sass_clear_options (struct Sass_Options* options);
  static void sass_reset_options (struct Sass_Options* options);
  static void copy_options(struct Sass_Options* to, struct Sass_Options* from) {
    // do not overwrite ourself
    if (to == from) return;
    // free assigned memory
    sass_clear_options(to);
    // move memory
    *to = *from;
    // Reset pointers on source
    sass_reset_options(from);
  }

  #define IMPLEMENT_SASS_OPTION_ACCESSOR(type, option) \
    type ADDCALL sass_option_get_##option (struct Sass_Options* options) { return options->option; } \
    void ADDCALL sass_option_set_##option (struct Sass_Options* options, type option) { options->option = option; }
  #define IMPLEMENT_SASS_OPTION_STRING_GETTER(type, option, def) \
    type ADDCALL sass_option_get_##option (struct Sass_Options* options) { return safe_str(options->option, def); }
  #define IMPLEMENT_SASS_OPTION_STRING_SETTER(type, option, def) \
    void ADDCALL sass_option_set_##option (struct Sass_Options* options, type option) \
    { free(options->option); options->option = option || def ? sass_copy_c_string(option ? option : def) : 0; }
  #define IMPLEMENT_SASS_OPTION_STRING_ACCESSOR(type, option, def) \
    IMPLEMENT_SASS_OPTION_STRING_GETTER(type, option, def) \
    IMPLEMENT_SASS_OPTION_STRING_SETTER(type, option, def)

  #define IMPLEMENT_SASS_CONTEXT_GETTER(type, option) \
    type ADDCALL sass_context_get_##option (struct Sass_Context* ctx) { return ctx->option; }
  #define IMPLEMENT_SASS_CONTEXT_TAKER(type, option) \
    type sass_context_take_##option (struct Sass_Context* ctx) \
    { type foo = ctx->option; ctx->option = 0; return foo; }


  // generic compilation function (not exported, use file/data compile instead)
  static Sass_Compiler* sass_prepare_context (Sass_Context* c_ctx, Context* cpp_ctx) throw()
  {
    try {
      // register our custom functions
      if (c_ctx->c_functions) {
        auto this_func_data = c_ctx->c_functions;
        while (this_func_data && *this_func_data) {
          cpp_ctx->add_c_function(*this_func_data);
          ++this_func_data;
        }
      }

      // register our custom headers
      if (c_ctx->c_headers) {
        auto this_head_data = c_ctx->c_headers;
        while (this_head_data && *this_head_data) {
          cpp_ctx->add_c_header(*this_head_data);
          ++this_head_data;
        }
      }

      // register our custom importers
      if (c_ctx->c_importers) {
        auto this_imp_data = c_ctx->c_importers;
        while (this_imp_data && *this_imp_data) {
          cpp_ctx->add_c_importer(*this_imp_data);
          ++this_imp_data;
        }
      }

      // reset error status
      c_ctx->error_json = 0;
      c_ctx->error_text = 0;
      c_ctx->error_message = 0;
      c_ctx->error_status = 0;
      // reset error position
      c_ctx->error_src = 0;
      c_ctx->error_file = 0;
      c_ctx->error_line = std::string::npos;
      c_ctx->error_column = std::string::npos;

      // allocate a new compiler instance
      void* ctxmem = calloc(1, sizeof(struct Sass_Compiler));
      if (ctxmem == 0) { std::cerr << "Error allocating memory for context" << std::endl; return 0; }
      Sass_Compiler* compiler = (struct Sass_Compiler*) ctxmem;
      compiler->state = SASS_COMPILER_CREATED;

      // store in sass compiler
      compiler->c_ctx = c_ctx;
      compiler->cpp_ctx = cpp_ctx;
      cpp_ctx->c_compiler = compiler;

      // use to parse block
      return compiler;

    }
    // pass errors to generic error handler
    catch (...) { handle_errors(c_ctx); }

    // error
    return 0;

  }

  // generic compilation function (not exported, use file/data compile instead)
  static int sass_compile_context (Sass_Context* c_ctx, Context* cpp_ctx)
  {

    // prepare sass compiler with context and options
    Sass_Compiler* compiler = sass_prepare_context(c_ctx, cpp_ctx);

    try {
      // call each compiler step
      sass_compiler_parse(compiler);
      sass_compiler_execute(compiler);
    }
    // pass errors to generic error handler
    catch (...) { handle_errors(c_ctx); }

    sass_delete_compiler(compiler);

    return c_ctx->error_status;
  }

  inline void init_options (struct Sass_Options* options)
  {
    options->precision = 10;
    options->indent = "  ";
    options->linefeed = LFEED;
  }

  Sass_Options* ADDCALL sass_make_options (void)
  {
    struct Sass_Options* options = (struct Sass_Options*) calloc(1, sizeof(struct Sass_Options));
    if (options == 0) { std::cerr << "Error allocating memory for options" << std::endl; return 0; }
    init_options(options);
    return options;
  }

  Sass_File_Context* ADDCALL sass_make_file_context(const char* input_path)
  {
    SharedObj::setTaint(true); // needed for static colors
    struct Sass_File_Context* ctx = (struct Sass_File_Context*) calloc(1, sizeof(struct Sass_File_Context));
    if (ctx == 0) { std::cerr << "Error allocating memory for file context" << std::endl; return 0; }
    ctx->type = SASS_CONTEXT_FILE;
    init_options(ctx);
    try {
      if (input_path == 0) { throw(std::runtime_error("File context created without an input path")); }
      if (*input_path == 0) { throw(std::runtime_error("File context created with empty input path")); }
      sass_option_set_input_path(ctx, input_path);
    } catch (...) {
      handle_errors(ctx);
    }
    return ctx;
  }

  Sass_Data_Context* ADDCALL sass_make_data_context(char* source_string)
  {
    struct Sass_Data_Context* ctx = (struct Sass_Data_Context*) calloc(1, sizeof(struct Sass_Data_Context));
    if (ctx == 0) { std::cerr << "Error allocating memory for data context" << std::endl; return 0; }
    ctx->type = SASS_CONTEXT_DATA;
    init_options(ctx);
    try {
      if (source_string == 0) { throw(std::runtime_error("Data context created without a source string")); }
      if (*source_string == 0) { throw(std::runtime_error("Data context created with empty source string")); }
      ctx->source_string = source_string;
    } catch (...) {
      handle_errors(ctx);
    }
    return ctx;
  }

  struct Sass_Compiler* ADDCALL sass_make_data_compiler (struct Sass_Data_Context* data_ctx)
  {
    if (data_ctx == 0) return 0;
    Context* cpp_ctx = new Data_Context(*data_ctx);
    return sass_prepare_context(data_ctx, cpp_ctx);
  }

  struct Sass_Compiler* ADDCALL sass_make_file_compiler (struct Sass_File_Context* file_ctx)
  {
    if (file_ctx == 0) return 0;
    Context* cpp_ctx = new File_Context(*file_ctx);
    return sass_prepare_context(file_ctx, cpp_ctx);
  }

  int ADDCALL sass_compile_data_context(Sass_Data_Context* data_ctx)
  {
    if (data_ctx == 0) return 1;
    if (data_ctx->error_status)
      return data_ctx->error_status;
    try {
      if (data_ctx->source_string == 0) { throw(std::runtime_error("Data context has no source string")); }
      // empty source string is a valid case, even if not really usefull (different than with file context)
      // if (*data_ctx->source_string == 0) { throw(std::runtime_error("Data context has empty source string")); }
    }
    catch (...) { return handle_errors(data_ctx) | 1; }
    Context* cpp_ctx = new Data_Context(*data_ctx);
    return sass_compile_context(data_ctx, cpp_ctx);
  }

  int ADDCALL sass_compile_file_context(Sass_File_Context* file_ctx)
  {
    if (file_ctx == 0) return 1;
    if (file_ctx->error_status)
      return file_ctx->error_status;
    try {
      if (file_ctx->input_path == 0) { throw(std::runtime_error("File context has no input path")); }
      if (*file_ctx->input_path == 0) { throw(std::runtime_error("File context has empty input path")); }
    }
    catch (...) { return handle_errors(file_ctx) | 1; }
    Context* cpp_ctx = new File_Context(*file_ctx);
    return sass_compile_context(file_ctx, cpp_ctx);
  }

  int ADDCALL sass_compiler_parse(struct Sass_Compiler* compiler)
  {
    if (compiler == 0) return 1;
    if (compiler->state == SASS_COMPILER_PARSED) return 0;
    if (compiler->state != SASS_COMPILER_CREATED) return -1;
    if (compiler->c_ctx == NULL) return 1;
    if (compiler->cpp_ctx == NULL) return 1;
    if (compiler->c_ctx->error_status)
      return compiler->c_ctx->error_status;
    // parse the context we have set up (file or data)
    compiler->root = sass_parse_block(compiler);
    // success
    return 0;
  }

  int ADDCALL sass_compiler_execute(struct Sass_Compiler* compiler)
  {
    if (compiler == 0) return 1;
    if (compiler->state == SASS_COMPILER_EXECUTED) return 0;
    if (compiler->state != SASS_COMPILER_PARSED) return -1;
    if (compiler->c_ctx == NULL) return 1;
    if (compiler->cpp_ctx == NULL) return 1;
    if (compiler->root.isNull()) return 1;
    if (compiler->c_ctx->error_status)
      return compiler->c_ctx->error_status;
    compiler->state = SASS_COMPILER_EXECUTED;
    Context* cpp_ctx = compiler->cpp_ctx;
    Block_Obj root = compiler->root;
    // compile the parsed root block
    try { compiler->c_ctx->output_string = cpp_ctx->render(root); }
    // pass catched errors to generic error handler
    catch (...) { return handle_errors(compiler->c_ctx) | 1; }
    // generate source map json and store on context
    compiler->c_ctx->source_map_string = cpp_ctx->render_srcmap();
    // success
    return 0;
  }

  // helper function, not exported, only accessible locally
  static void sass_reset_options (struct Sass_Options* options)
  {
    // free pointer before
    // or copy/move them
    options->input_path = 0;
    options->output_path = 0;
    options->plugin_path = 0;
    options->include_path = 0;
    options->source_map_file = 0;
    options->source_map_root = 0;
    options->c_functions = 0;
    options->c_importers = 0;
    options->c_headers = 0;
    options->plugin_paths = 0;
    options->include_paths = 0;
  }

  // helper function, not exported, only accessible locally
  static void sass_clear_options (struct Sass_Options* options)
  {
    if (options == 0) return;
    // Deallocate custom functions, headers and importes
    sass_delete_function_list(options->c_functions);
    sass_delete_importer_list(options->c_importers);
    sass_delete_importer_list(options->c_headers);
    // Deallocate inc paths
    if (options->plugin_paths) {
      struct string_list* cur;
      struct string_list* next;
      cur = options->plugin_paths;
      while (cur) {
        next = cur->next;
        free(cur->string);
        free(cur);
        cur = next;
      }
    }
    // Deallocate inc paths
    if (options->include_paths) {
      struct string_list* cur;
      struct string_list* next;
      cur = options->include_paths;
      while (cur) {
        next = cur->next;
        free(cur->string);
        free(cur);
        cur = next;
      }
    }
    // Free options strings
    free(options->input_path);
    free(options->output_path);
    free(options->plugin_path);
    free(options->include_path);
    free(options->source_map_file);
    free(options->source_map_root);
    // Reset our pointers
    options->input_path = 0;
    options->output_path = 0;
    options->plugin_path = 0;
    options->include_path = 0;
    options->source_map_file = 0;
    options->source_map_root = 0;
    options->c_functions = 0;
    options->c_importers = 0;
    options->c_headers = 0;
    options->plugin_paths = 0;
    options->include_paths = 0;
  }

  // helper function, not exported, only accessible locally
  // sass_free_context is also defined in old sass_interface
  static void sass_clear_context (struct Sass_Context* ctx)
  {
    if (ctx == 0) return;
    // release the allocated memory (mostly via sass_copy_c_string)
    if (ctx->output_string)     free(ctx->output_string);
    if (ctx->source_map_string) free(ctx->source_map_string);
    if (ctx->error_message)     free(ctx->error_message);
    if (ctx->error_text)        free(ctx->error_text);
    if (ctx->error_json)        free(ctx->error_json);
    if (ctx->error_file)        free(ctx->error_file);
    free_string_array(ctx->included_files);
    // play safe and reset properties
    ctx->output_string = 0;
    ctx->source_map_string = 0;
    ctx->error_message = 0;
    ctx->error_text = 0;
    ctx->error_json = 0;
    ctx->error_file = 0;
    ctx->included_files = 0;
    // debug leaked memory
    #ifdef DEBUG_SHARED_PTR
      SharedObj::dumpMemLeaks();
    #endif
    // now clear the options
    sass_clear_options(ctx);
  }

  void ADDCALL sass_delete_compiler (struct Sass_Compiler* compiler)
  {
    if (compiler == 0) {
      return;
    }
    Context* cpp_ctx = compiler->cpp_ctx;
    if (cpp_ctx) delete(cpp_ctx);
    compiler->cpp_ctx = NULL;
    compiler->c_ctx = NULL;
    compiler->root = {};
    free(compiler);
  }

  void ADDCALL sass_delete_options (struct Sass_Options* options)
  {
    sass_clear_options(options); free(options);
  }

  // Deallocate all associated memory with file context
  void ADDCALL sass_delete_file_context (struct Sass_File_Context* ctx)
  {
    // clear the context and free it
    sass_clear_context(ctx); free(ctx);
  }
  // Deallocate all associated memory with data context
  void ADDCALL sass_delete_data_context (struct Sass_Data_Context* ctx)
  {
    // clean the source string if it was not passed
    // we reset this member once we start parsing
    if (ctx->source_string) free(ctx->source_string);
    // clear the context and free it
    sass_clear_context(ctx); free(ctx);
  }

  // Getters for sass context from specific implementations
  struct Sass_Context* ADDCALL sass_file_context_get_context(struct Sass_File_Context* ctx) { return ctx; }
  struct Sass_Context* ADDCALL sass_data_context_get_context(struct Sass_Data_Context* ctx) { return ctx; }

  // Getters for context options from Sass_Context
  struct Sass_Options* ADDCALL sass_context_get_options(struct Sass_Context* ctx) { return ctx; }
  struct Sass_Options* ADDCALL sass_file_context_get_options(struct Sass_File_Context* ctx) { return ctx; }
  struct Sass_Options* ADDCALL sass_data_context_get_options(struct Sass_Data_Context* ctx) { return ctx; }
  void ADDCALL sass_file_context_set_options (struct Sass_File_Context* ctx, struct Sass_Options* opt) { copy_options(ctx, opt); }
  void ADDCALL sass_data_context_set_options (struct Sass_Data_Context* ctx, struct Sass_Options* opt) { copy_options(ctx, opt); }

  // Getters for Sass_Compiler options (get conected sass context)
  enum Sass_Compiler_State ADDCALL sass_compiler_get_state(struct Sass_Compiler* compiler) { return compiler->state; }
  struct Sass_Context* ADDCALL sass_compiler_get_context(struct Sass_Compiler* compiler) { return compiler->c_ctx; }
  struct Sass_Options* ADDCALL sass_compiler_get_options(struct Sass_Compiler* compiler) { return compiler->c_ctx; }
  // Getters for Sass_Compiler options (query import stack)
  size_t ADDCALL sass_compiler_get_import_stack_size(struct Sass_Compiler* compiler) { return compiler->cpp_ctx->import_stack.size(); }
  Sass_Import_Entry ADDCALL sass_compiler_get_last_import(struct Sass_Compiler* compiler) { return compiler->cpp_ctx->import_stack.back(); }
  Sass_Import_Entry ADDCALL sass_compiler_get_import_entry(struct Sass_Compiler* compiler, size_t idx) { return compiler->cpp_ctx->import_stack[idx]; }
  // Getters for Sass_Compiler options (query function stack)
  size_t ADDCALL sass_compiler_get_callee_stack_size(struct Sass_Compiler* compiler) { return compiler->cpp_ctx->callee_stack.size(); }
  Sass_Callee_Entry ADDCALL sass_compiler_get_last_callee(struct Sass_Compiler* compiler) { return &compiler->cpp_ctx->callee_stack.back(); }
  Sass_Callee_Entry ADDCALL sass_compiler_get_callee_entry(struct Sass_Compiler* compiler, size_t idx) { return &compiler->cpp_ctx->callee_stack[idx]; }

  // Calculate the size of the stored null terminated array
  size_t ADDCALL sass_context_get_included_files_size (struct Sass_Context* ctx)
  { size_t l = 0; auto i = ctx->included_files; while (i && *i) { ++i; ++l; } return l; }

  // Create getter and setters for options
  IMPLEMENT_SASS_OPTION_ACCESSOR(int, precision);
  IMPLEMENT_SASS_OPTION_ACCESSOR(enum Sass_Output_Style, output_style);
  IMPLEMENT_SASS_OPTION_ACCESSOR(bool, source_comments);
  IMPLEMENT_SASS_OPTION_ACCESSOR(bool, source_map_embed);
  IMPLEMENT_SASS_OPTION_ACCESSOR(bool, source_map_contents);
  IMPLEMENT_SASS_OPTION_ACCESSOR(bool, source_map_file_urls);
  IMPLEMENT_SASS_OPTION_ACCESSOR(bool, omit_source_map_url);
  IMPLEMENT_SASS_OPTION_ACCESSOR(bool, is_indented_syntax_src);
  IMPLEMENT_SASS_OPTION_ACCESSOR(Sass_Function_List, c_functions);
  IMPLEMENT_SASS_OPTION_ACCESSOR(Sass_Importer_List, c_importers);
  IMPLEMENT_SASS_OPTION_ACCESSOR(Sass_Importer_List, c_headers);
  IMPLEMENT_SASS_OPTION_ACCESSOR(const char*, indent);
  IMPLEMENT_SASS_OPTION_ACCESSOR(const char*, linefeed);
  IMPLEMENT_SASS_OPTION_STRING_SETTER(const char*, plugin_path, 0);
  IMPLEMENT_SASS_OPTION_STRING_SETTER(const char*, include_path, 0);
  IMPLEMENT_SASS_OPTION_STRING_ACCESSOR(const char*, input_path, 0);
  IMPLEMENT_SASS_OPTION_STRING_ACCESSOR(const char*, output_path, 0);
  IMPLEMENT_SASS_OPTION_STRING_ACCESSOR(const char*, source_map_file, 0);
  IMPLEMENT_SASS_OPTION_STRING_ACCESSOR(const char*, source_map_root, 0);

  // Create getter and setters for context
  IMPLEMENT_SASS_CONTEXT_GETTER(int, error_status);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, error_json);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, error_message);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, error_text);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, error_file);
  IMPLEMENT_SASS_CONTEXT_GETTER(size_t, error_line);
  IMPLEMENT_SASS_CONTEXT_GETTER(size_t, error_column);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, error_src);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, output_string);
  IMPLEMENT_SASS_CONTEXT_GETTER(const char*, source_map_string);
  IMPLEMENT_SASS_CONTEXT_GETTER(char**, included_files);

  // Take ownership of memory (value on context is set to 0)
  IMPLEMENT_SASS_CONTEXT_TAKER(char*, error_json);
  IMPLEMENT_SASS_CONTEXT_TAKER(char*, error_message);
  IMPLEMENT_SASS_CONTEXT_TAKER(char*, error_text);
  IMPLEMENT_SASS_CONTEXT_TAKER(char*, error_file);
  IMPLEMENT_SASS_CONTEXT_TAKER(char*, output_string);
  IMPLEMENT_SASS_CONTEXT_TAKER(char*, source_map_string);
  IMPLEMENT_SASS_CONTEXT_TAKER(char**, included_files);

  // Push function for include paths (no manipulation support for now)
  void ADDCALL sass_option_push_include_path(struct Sass_Options* options, const char* path)
  {

    struct string_list* include_path = (struct string_list*) calloc(1, sizeof(struct string_list));
    if (include_path == 0) return;
    include_path->string = path ? sass_copy_c_string(path) : 0;
    struct string_list* last = options->include_paths;
    if (!options->include_paths) {
      options->include_paths = include_path;
    } else {
      while (last->next)
        last = last->next;
      last->next = include_path;
    }

  }

  // Push function for include paths (no manipulation support for now)
  size_t ADDCALL sass_option_get_include_path_size(struct Sass_Options* options)
  {
    size_t len = 0;
    struct string_list* cur = options->include_paths;
    while (cur) { len ++; cur = cur->next; }
    return len;
  }

  // Push function for include paths (no manipulation support for now)
  const char* ADDCALL sass_option_get_include_path(struct Sass_Options* options, size_t i)
  {
    struct string_list* cur = options->include_paths;
    while (i) { i--; cur = cur->next; }
    return cur->string;
  }

  // Push function for plugin paths (no manipulation support for now)
  void ADDCALL sass_option_push_plugin_path(struct Sass_Options* options, const char* path)
  {

    struct string_list* plugin_path = (struct string_list*) calloc(1, sizeof(struct string_list));
    if (plugin_path == 0) return;
    plugin_path->string = path ? sass_copy_c_string(path) : 0;
    struct string_list* last = options->plugin_paths;
    if (!options->plugin_paths) {
      options->plugin_paths = plugin_path;
    } else {
      while (last->next)
        last = last->next;
      last->next = plugin_path;
    }

  }

}

namespace Sass {

  #define IMPLEMENT_BASE_CAST(T) \
  template<> \
  T* Cast(AST_Node* ptr) { \
    return dynamic_cast<T*>(ptr); \
  }; \
  \
  template<> \
  const T* Cast(const AST_Node* ptr) { \
    return dynamic_cast<const T*>(ptr); \
  }; \

  IMPLEMENT_BASE_CAST(AST_Node)
  IMPLEMENT_BASE_CAST(Expression)
  IMPLEMENT_BASE_CAST(Statement)
  IMPLEMENT_BASE_CAST(Has_Block)
  IMPLEMENT_BASE_CAST(PreValue)
  IMPLEMENT_BASE_CAST(Value)
  IMPLEMENT_BASE_CAST(Color)
  IMPLEMENT_BASE_CAST(List)
  IMPLEMENT_BASE_CAST(String)
  IMPLEMENT_BASE_CAST(String_Constant)
  IMPLEMENT_BASE_CAST(Supports_Condition)
  IMPLEMENT_BASE_CAST(Selector)
  IMPLEMENT_BASE_CAST(Simple_Selector)

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  namespace Exception {

    Base::Base(ParserState pstate, std::string msg, Backtraces traces)
    : std::runtime_error(msg), msg(msg),
      prefix("Error"), pstate(pstate), traces(traces)
    { }

    InvalidSass::InvalidSass(ParserState pstate, Backtraces traces, std::string msg, char* owned_src)
    : Base(pstate, msg, traces), owned_src(owned_src)
    { }


    InvalidParent::InvalidParent(Selector_Ptr parent, Backtraces traces, Selector_Ptr selector)
    : Base(selector->pstate(), def_msg, traces), parent(parent), selector(selector)
    {
      msg = "Invalid parent selector for "
        "\"" + selector->to_string(Sass_Inspect_Options()) + "\": "
        "\"" + parent->to_string(Sass_Inspect_Options()) + "\"";
    }

    InvalidVarKwdType::InvalidVarKwdType(ParserState pstate, Backtraces traces, std::string name, const Argument_Ptr arg)
    : Base(pstate, def_msg, traces), name(name), arg(arg)
    {
      msg = "Variable keyword argument map must have string keys.\n" +
        name + " is not a string in " + arg->to_string() + ".";
    }

    InvalidArgumentType::InvalidArgumentType(ParserState pstate, Backtraces traces, std::string fn, std::string arg, std::string type, const Value_Ptr value)
    : Base(pstate, def_msg, traces), fn(fn), arg(arg), type(type), value(value)
    {
      msg = arg + ": \"";
      if (value) msg += value->to_string(Sass_Inspect_Options());
      msg += "\" is not a " + type + " for `" + fn + "'";
    }

    MissingArgument::MissingArgument(ParserState pstate, Backtraces traces, std::string fn, std::string arg, std::string fntype)
    : Base(pstate, def_msg, traces), fn(fn), arg(arg), fntype(fntype)
    {
      msg = fntype + " " + fn + " is missing argument " + arg + ".";
    }

    InvalidSyntax::InvalidSyntax(ParserState pstate, Backtraces traces, std::string msg)
    : Base(pstate, msg, traces)
    { }

    NestingLimitError::NestingLimitError(ParserState pstate, Backtraces traces, std::string msg)
    : Base(pstate, msg, traces)
    { }

    DuplicateKeyError::DuplicateKeyError(Backtraces traces, const Map& dup, const Expression& org)
    : Base(org.pstate(), def_msg, traces), dup(dup), org(org)
    {
      msg = "Duplicate key " + dup.get_duplicate_key()->inspect() + " in map (" + org.inspect() + ").";
    }

    TypeMismatch::TypeMismatch(Backtraces traces, const Expression& var, const std::string type)
    : Base(var.pstate(), def_msg, traces), var(var), type(type)
    {
      msg = var.to_string() + " is not an " + type + ".";
    }

    InvalidValue::InvalidValue(Backtraces traces, const Expression& val)
    : Base(val.pstate(), def_msg, traces), val(val)
    {
      msg = val.to_string() + " isn't a valid CSS value.";
    }

    StackError::StackError(Backtraces traces, const AST_Node& node)
    : Base(node.pstate(), def_msg, traces), node(node)
    {
      msg = "stack level too deep";
    }

    IncompatibleUnits::IncompatibleUnits(const Units& lhs, const Units& rhs)
    {
      msg = "Incompatible units: '" + rhs.unit() + "' and '" + lhs.unit() + "'.";
    }

    IncompatibleUnits::IncompatibleUnits(const UnitType lhs, const UnitType rhs)
    {
      msg = std::string("Incompatible units: '") + unit_to_string(rhs) + "' and '" + unit_to_string(lhs) + "'.";
    }

    AlphaChannelsNotEqual::AlphaChannelsNotEqual(Expression_Ptr_Const lhs, Expression_Ptr_Const rhs, enum Sass_OP op)
    : OperationError(), lhs(lhs), rhs(rhs), op(op)
    {
      msg = "Alpha channels must be equal: " +
        lhs->to_string({ NESTED, 5 }) +
        " " + sass_op_to_name(op) + " " +
        rhs->to_string({ NESTED, 5 }) + ".";
    }

    ZeroDivisionError::ZeroDivisionError(const Expression& lhs, const Expression& rhs)
    : OperationError(), lhs(lhs), rhs(rhs)
    {
      msg = "divided by 0";
    }

    UndefinedOperation::UndefinedOperation(Expression_Ptr_Const lhs, Expression_Ptr_Const rhs, enum Sass_OP op)
    : OperationError(), lhs(lhs), rhs(rhs), op(op)
    {
      msg = def_op_msg + ": \"" +
        lhs->to_string({ NESTED, 5 }) +
        " " + sass_op_to_name(op) + " " +
        rhs->to_string({ TO_SASS, 5 }) +
        "\".";
    }

    InvalidNullOperation::InvalidNullOperation(Expression_Ptr_Const lhs, Expression_Ptr_Const rhs, enum Sass_OP op)
    : UndefinedOperation(lhs, rhs, op)
    {
      msg = def_op_null_msg + ": \"" + lhs->inspect() + " " + sass_op_to_name(op) + " " + rhs->inspect() + "\".";
    }

    SassValueError::SassValueError(Backtraces traces, ParserState pstate, OperationError& err)
    : Base(pstate, err.what(), traces)
    {
      msg = err.what();
      prefix = err.errtype();
    }

  }


  void warn(std::string msg, ParserState pstate)
  {
    std::cerr << "Warning: " << msg << std::endl;
  }

  void warning(std::string msg, ParserState pstate)
  {
    std::string cwd(Sass::File::get_cwd());
    std::string abs_path(Sass::File::rel2abs(pstate.path, cwd, cwd));
    std::string rel_path(Sass::File::abs2rel(pstate.path, cwd, cwd));
    std::string output_path(Sass::File::path_for_console(rel_path, abs_path, pstate.path));

    std::cerr << "WARNING on line " << pstate.line+1 << ", column " << pstate.column+1 << " of " << output_path << ":" << std::endl;
    std::cerr << msg << std::endl << std::endl;
  }

  void warn(std::string msg, ParserState pstate, Backtrace* bt)
  {
    warn(msg, pstate);
  }

  void deprecated_function(std::string msg, ParserState pstate)
  {
    std::string cwd(Sass::File::get_cwd());
    std::string abs_path(Sass::File::rel2abs(pstate.path, cwd, cwd));
    std::string rel_path(Sass::File::abs2rel(pstate.path, cwd, cwd));
    std::string output_path(Sass::File::path_for_console(rel_path, abs_path, pstate.path));

    std::cerr << "DEPRECATION WARNING: " << msg << std::endl;
    std::cerr << "will be an error in future versions of Sass." << std::endl;
    std::cerr << "        on line " << pstate.line+1 << " of " << output_path << std::endl;
  }

  void deprecated(std::string msg, std::string msg2, bool with_column, ParserState pstate)
  {
    std::string cwd(Sass::File::get_cwd());
    std::string abs_path(Sass::File::rel2abs(pstate.path, cwd, cwd));
    std::string rel_path(Sass::File::abs2rel(pstate.path, cwd, cwd));
    std::string output_path(Sass::File::path_for_console(rel_path, pstate.path, pstate.path));

    std::cerr << "DEPRECATION WARNING on line " << pstate.line + 1;
    if (with_column) std::cerr << ", column " << pstate.column + pstate.offset.column + 1;
    if (output_path.length()) std::cerr << " of " << output_path;
    std::cerr << ":" << std::endl;
    std::cerr << msg << std::endl;
    if (msg2.length()) std::cerr << msg2 << std::endl;
    std::cerr << std::endl;
  }

  void deprecated_bind(std::string msg, ParserState pstate)
  {
    std::string cwd(Sass::File::get_cwd());
    std::string abs_path(Sass::File::rel2abs(pstate.path, cwd, cwd));
    std::string rel_path(Sass::File::abs2rel(pstate.path, cwd, cwd));
    std::string output_path(Sass::File::path_for_console(rel_path, abs_path, pstate.path));

    std::cerr << "WARNING: " << msg << std::endl;
    std::cerr << "        on line " << pstate.line+1 << " of " << output_path << std::endl;
    std::cerr << "This will be an error in future versions of Sass." << std::endl;
  }

  // should be replaced with error with backtraces
  void coreError(std::string msg, ParserState pstate)
  {
    Backtraces traces;
    throw Exception::InvalidSyntax(pstate, traces, msg);
  }

  void error(std::string msg, ParserState pstate, Backtraces& traces)
  {
    traces.push_back(Backtrace(pstate));
    throw Exception::InvalidSyntax(pstate, traces, msg);
  }

}
#include <array>

namespace Sass {

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Selector::Selector(ParserState pstate)
  : Expression(pstate),
    has_line_feed_(false),
    has_line_break_(false),
    is_optional_(false),
    media_block_(0),
    hash_(0)
  { concrete_type(SELECTOR); }

  Selector::Selector(const Selector* ptr)
  : Expression(ptr),
    has_line_feed_(ptr->has_line_feed_),
    has_line_break_(ptr->has_line_break_),
    is_optional_(ptr->is_optional_),
    media_block_(ptr->media_block_),
    hash_(ptr->hash_)
  { concrete_type(SELECTOR); }

  void Selector::set_media_block(Media_Block_Ptr mb)
  {
    media_block(mb);
  }

  bool Selector::has_parent_ref() const
  {
    return false;
  }

  bool Selector::has_real_parent_ref() const
  {
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Selector_Schema::Selector_Schema(ParserState pstate, String_Obj c)
  : AST_Node(pstate),
    contents_(c),
    connect_parent_(true),
    media_block_(NULL),
    hash_(0)
  { }
  Selector_Schema::Selector_Schema(const Selector_Schema* ptr)
  : AST_Node(ptr),
    contents_(ptr->contents_),
    connect_parent_(ptr->connect_parent_),
    media_block_(ptr->media_block_),
    hash_(ptr->hash_)
  { }

  unsigned long Selector_Schema::specificity() const
  {
    return 0;
  }

  size_t Selector_Schema::hash() const {
    if (hash_ == 0) {
      hash_combine(hash_, contents_->hash());
    }
    return hash_;
  }

  bool Selector_Schema::has_parent_ref() const
  {
    if (String_Schema_Obj schema = Cast<String_Schema>(contents())) {
      if (schema->empty()) return false;
      const auto& first = *schema->at(0);
      return typeid(first) == typeid(Parent_Selector);
    }
    return false;
  }

  bool Selector_Schema::has_real_parent_ref() const
  {
    if (String_Schema_Obj schema = Cast<String_Schema>(contents())) {
      if (schema->empty()) return false;
      const auto& first = *schema->at(0);
      return typeid(first) == typeid(Parent_Reference);
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Simple_Selector::Simple_Selector(ParserState pstate, std::string n)
  : Selector(pstate), ns_(""), name_(n), has_ns_(false)
  {
    size_t pos = n.find('|');
    // found some namespace
    if (pos != std::string::npos) {
      has_ns_ = true;
      ns_ = n.substr(0, pos);
      name_ = n.substr(pos + 1);
    }
  }
  Simple_Selector::Simple_Selector(const Simple_Selector* ptr)
  : Selector(ptr),
    ns_(ptr->ns_),
    name_(ptr->name_),
    has_ns_(ptr->has_ns_)
  { }

  std::string Simple_Selector::ns_name() const
  {
    std::string name("");
    if (has_ns_)
      name += ns_ + "|";
    return name + name_;
  }

  size_t Simple_Selector::hash() const
  {
    if (hash_ == 0) {
      hash_combine(hash_, std::hash<int>()(SELECTOR));
      hash_combine(hash_, std::hash<int>()(simple_type()));
      if (!name_.empty()) hash_combine(hash_, std::hash<std::string>()(name()));
      if (has_ns_) hash_combine(hash_, std::hash<std::string>()(ns()));
    }
    return hash_;
  }

  bool Simple_Selector::empty() const {
    return ns().empty() && name().empty();
  }

  // namespace compare functions
  bool Simple_Selector::is_ns_eq(const Simple_Selector& r) const
  {
    return has_ns_ == r.has_ns_ && ns_ == r.ns_;
  }

  // namespace query functions
  bool Simple_Selector::is_universal_ns() const
  {
    return has_ns_ && ns_ == "*";
  }

  bool Simple_Selector::is_empty_ns() const
  {
    return !has_ns_ || ns_ == "";
  }

  bool Simple_Selector::has_empty_ns() const
  {
    return has_ns_ && ns_ == "";
  }

  bool Simple_Selector::has_qualified_ns() const
  {
    return has_ns_ && ns_ != "" && ns_ != "*";
  }

  // name query functions
  bool Simple_Selector::is_universal() const
  {
    return name_ == "*";
  }

  bool Simple_Selector::has_placeholder()
  {
    return false;
  }

  bool Simple_Selector::has_parent_ref() const
  {
    return false;
  };

  bool Simple_Selector::has_real_parent_ref() const
  {
    return false;
  };

  bool Simple_Selector::is_pseudo_element() const
  {
    return false;
  }

  bool Simple_Selector::is_superselector_of(Compound_Selector_Ptr_Const sub) const
  {
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Parent_Selector::Parent_Selector(ParserState pstate, bool r)
  : Simple_Selector(pstate, "&"), real_(r)
  { simple_type(PARENT_SEL); }
  Parent_Selector::Parent_Selector(const Parent_Selector* ptr)
  : Simple_Selector(ptr), real_(ptr->real_)
  { simple_type(PARENT_SEL); }

  bool Parent_Selector::has_parent_ref() const
  {
    return true;
  };

  bool Parent_Selector::has_real_parent_ref() const
  {
    return real();
  };

  unsigned long Parent_Selector::specificity() const
  {
    return 0;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Placeholder_Selector::Placeholder_Selector(ParserState pstate, std::string n)
  : Simple_Selector(pstate, n)
  { simple_type(PLACEHOLDER_SEL); }
  Placeholder_Selector::Placeholder_Selector(const Placeholder_Selector* ptr)
  : Simple_Selector(ptr)
  { simple_type(PLACEHOLDER_SEL); }
  unsigned long Placeholder_Selector::specificity() const
  {
    return Constants::Specificity_Base;
  }
  bool Placeholder_Selector::has_placeholder() {
    return true;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Type_Selector::Type_Selector(ParserState pstate, std::string n)
  : Simple_Selector(pstate, n)
  { simple_type(TYPE_SEL); }
  Type_Selector::Type_Selector(const Type_Selector* ptr)
  : Simple_Selector(ptr)
  { simple_type(TYPE_SEL); }

  unsigned long Type_Selector::specificity() const
  {
    if (name() == "*") return 0;
    else return Constants::Specificity_Element;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Class_Selector::Class_Selector(ParserState pstate, std::string n)
  : Simple_Selector(pstate, n)
  { simple_type(CLASS_SEL); }
  Class_Selector::Class_Selector(const Class_Selector* ptr)
  : Simple_Selector(ptr)
  { simple_type(CLASS_SEL); }

  unsigned long Class_Selector::specificity() const
  {
    return Constants::Specificity_Class;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Id_Selector::Id_Selector(ParserState pstate, std::string n)
  : Simple_Selector(pstate, n)
  { simple_type(ID_SEL); }
  Id_Selector::Id_Selector(const Id_Selector* ptr)
  : Simple_Selector(ptr)
  { simple_type(ID_SEL); }

  unsigned long Id_Selector::specificity() const
  {
    return Constants::Specificity_ID;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Attribute_Selector::Attribute_Selector(ParserState pstate, std::string n, std::string m, String_Obj v, char o)
  : Simple_Selector(pstate, n), matcher_(m), value_(v), modifier_(o)
  { simple_type(ATTRIBUTE_SEL); }
  Attribute_Selector::Attribute_Selector(const Attribute_Selector* ptr)
  : Simple_Selector(ptr),
    matcher_(ptr->matcher_),
    value_(ptr->value_),
    modifier_(ptr->modifier_)
  { simple_type(ATTRIBUTE_SEL); }

  size_t Attribute_Selector::hash() const
  {
    if (hash_ == 0) {
      hash_combine(hash_, Simple_Selector::hash());
      hash_combine(hash_, std::hash<std::string>()(matcher()));
      if (value_) hash_combine(hash_, value_->hash());
    }
    return hash_;
  }

  unsigned long Attribute_Selector::specificity() const
  {
    return Constants::Specificity_Attr;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Pseudo_Selector::Pseudo_Selector(ParserState pstate, std::string n, String_Obj expr)
  : Simple_Selector(pstate, n), expression_(expr)
  { simple_type(PSEUDO_SEL); }
  Pseudo_Selector::Pseudo_Selector(const Pseudo_Selector* ptr)
  : Simple_Selector(ptr), expression_(ptr->expression_)
  { simple_type(PSEUDO_SEL); }

  // A pseudo-element is made of two colons (::) followed by the name.
  // The `::` notation is introduced by the current document in order to
  // establish a discrimination between pseudo-classes and pseudo-elements.
  // For compatibility with existing style sheets, user agents must also
  // accept the previous one-colon notation for pseudo-elements introduced
  // in CSS levels 1 and 2 (namely, :first-line, :first-letter, :before and
  // :after). This compatibility is not allowed for the new pseudo-elements
  // introduced in this specification.
  bool Pseudo_Selector::is_pseudo_element() const
  {
    return (name_[0] == ':' && name_[1] == ':')
            || is_pseudo_class_element(name_);
  }

  size_t Pseudo_Selector::hash() const
  {
    if (hash_ == 0) {
      hash_combine(hash_, Simple_Selector::hash());
      if (expression_) hash_combine(hash_, expression_->hash());
    }
    return hash_;
  }

  unsigned long Pseudo_Selector::specificity() const
  {
    if (is_pseudo_element())
      return Constants::Specificity_Element;
    return Constants::Specificity_Pseudo;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Wrapped_Selector::Wrapped_Selector(ParserState pstate, std::string n, Selector_List_Obj sel)
  : Simple_Selector(pstate, n), selector_(sel)
  { simple_type(WRAPPED_SEL); }
  Wrapped_Selector::Wrapped_Selector(const Wrapped_Selector* ptr)
  : Simple_Selector(ptr), selector_(ptr->selector_)
  { simple_type(WRAPPED_SEL); }

  bool Wrapped_Selector::is_superselector_of(Wrapped_Selector_Ptr_Const sub) const
  {
    if (this->name() != sub->name()) return false;
    if (this->name() == ":current") return false;
    if (Selector_List_Obj rhs_list = Cast<Selector_List>(sub->selector())) {
      if (Selector_List_Obj lhs_list = Cast<Selector_List>(selector())) {
        return lhs_list->is_superselector_of(rhs_list);
      }
    }
    coreError("is_superselector expected a Selector_List", sub->pstate());
    return false;
  }

  // Selectors inside the negation pseudo-class are counted like any
  // other, but the negation itself does not count as a pseudo-class.

  void Wrapped_Selector::cloneChildren()
  {
    selector(SASS_MEMORY_CLONE(selector()));
  }

  size_t Wrapped_Selector::hash() const
  {
    if (hash_ == 0) {
      hash_combine(hash_, Simple_Selector::hash());
      if (selector_) hash_combine(hash_, selector_->hash());
    }
    return hash_;
  }

  bool Wrapped_Selector::has_parent_ref() const {
    if (!selector()) return false;
    return selector()->has_parent_ref();
  }

  bool Wrapped_Selector::has_real_parent_ref() const {
    if (!selector()) return false;
    return selector()->has_real_parent_ref();
  }

  unsigned long Wrapped_Selector::specificity() const
  {
    return selector_ ? selector_->specificity() : 0;
  }

  bool Wrapped_Selector::find ( bool (*f)(AST_Node_Obj) )
  {
    // check children first
    if (selector_) {
      if (selector_->find(f)) return true;
    }
    // execute last
    return f(this);
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Compound_Selector::Compound_Selector(ParserState pstate, size_t s)
  : Selector(pstate),
    Vectorized<Simple_Selector_Obj>(s),
    extended_(false),
    has_parent_reference_(false)
  { }

  Compound_Selector::Compound_Selector(const Compound_Selector* ptr)
  : Selector(ptr),
    Vectorized<Simple_Selector_Obj>(*ptr),
    extended_(ptr->extended_),
    has_parent_reference_(ptr->has_parent_reference_)
  { }

  bool Compound_Selector::contains_placeholder() {
    for (size_t i = 0, L = length(); i < L; ++i) {
      if ((*this)[i]->has_placeholder()) return true;
    }
    return false;
  };

  void Compound_Selector::cloneChildren()
  {
    for (size_t i = 0, l = length(); i < l; i++) {
      at(i) = SASS_MEMORY_CLONE(at(i));
    }
  }

  bool Compound_Selector::find ( bool (*f)(AST_Node_Obj) )
  {
    // check children first
    for (Simple_Selector_Obj sel : elements()) {
      if (sel->find(f)) return true;
    }
    // execute last
    return f(this);
  }

  bool Compound_Selector::has_parent_ref() const
  {
    for (Simple_Selector_Obj s : *this) {
      if (s && s->has_parent_ref()) return true;
    }
    return false;
  }

  bool Compound_Selector::has_real_parent_ref() const
  {
    for (Simple_Selector_Obj s : *this) {
      if (s && s->has_real_parent_ref()) return true;
    }
    return false;
  }

  bool Compound_Selector::is_superselector_of(Selector_List_Ptr_Const rhs, std::string wrapped) const
  {
    for (Complex_Selector_Obj item : rhs->elements()) {
      if (is_superselector_of(item, wrapped)) return true;
    }
    return false;
  }

  bool Compound_Selector::is_superselector_of(Complex_Selector_Ptr_Const rhs, std::string wrapped) const
  {
    if (rhs->head()) return is_superselector_of(rhs->head(), wrapped);
    return false;
  }

  bool Compound_Selector::is_superselector_of(Compound_Selector_Ptr_Const rhs, std::string wrapping) const
  {
    // Check if pseudo-elements are the same between the selectors
    {
      std::array<std::set<std::string>, 2> pseudosets;
      std::array<Compound_Selector_Ptr_Const, 2> compounds = {{this, rhs}};
      for (int i = 0; i < 2; ++i) {
        for (const Simple_Selector_Obj& el : compounds[i]->elements()) {
          if (el->is_pseudo_element()) {
            std::string pseudo(el->to_string());
            // strip off colons to ensure :after matches ::after since ruby sass is forgiving
            pseudosets[i].insert(pseudo.substr(pseudo.find_first_not_of(":")));
          }
        }
      }
      if (pseudosets[0] != pseudosets[1]) return false;
    }

    {
      Simple_Selector_Ptr_Const lbase = this->base();
      Simple_Selector_Ptr_Const rbase = rhs->base();
      if (lbase && rbase) {
        return *lbase == *rbase &&
               contains_all(std::unordered_set<Simple_Selector_Ptr_Const, HashPtr, ComparePtrs>(rhs->begin(), rhs->end()),
                            std::unordered_set<Simple_Selector_Ptr_Const, HashPtr, ComparePtrs>(this->begin(), this->end()));
      }
    }

    std::unordered_set<Selector_Ptr_Const, HashPtr, ComparePtrs> lset;
    for (size_t i = 0, iL = length(); i < iL; ++i)
    {
      Selector_Ptr_Const wlhs = (*this)[i].ptr();
      // very special case for wrapped matches selector
      if (Wrapped_Selector_Ptr_Const wrapped = Cast<Wrapped_Selector>(wlhs)) {
        if (wrapped->name() == ":not") {
          if (Selector_List_Obj not_list = Cast<Selector_List>(wrapped->selector())) {
            if (not_list->is_superselector_of(rhs, wrapped->name())) return false;
          } else {
            throw std::runtime_error("wrapped not selector is not a list");
          }
        }
        if (wrapped->name() == ":matches" || wrapped->name() == ":-moz-any") {
          wlhs = wrapped->selector();
          if (Selector_List_Obj list = Cast<Selector_List>(wrapped->selector())) {
            if (Compound_Selector_Ptr_Const comp = Cast<Compound_Selector>(rhs)) {
              if (!wrapping.empty() && wrapping != wrapped->name()) return false;
              if (wrapping.empty() || wrapping != wrapped->name()) {;
                if (list->is_superselector_of(comp, wrapped->name())) return true;
              }
            }
          }
        }
        Simple_Selector_Ptr rhs_sel = nullptr;
        if (rhs->elements().size() > i) rhs_sel = (*rhs)[i];
        if (Wrapped_Selector_Ptr wrapped_r = Cast<Wrapped_Selector>(rhs_sel)) {
          if (wrapped->name() == wrapped_r->name()) {
          if (wrapped->is_superselector_of(wrapped_r)) {
             continue;
          }}
        }
      }
      lset.insert(wlhs);
    }

    if (lset.empty()) return true;

    std::unordered_set<Selector_Ptr_Const, HashPtr, ComparePtrs> rset;
    for (size_t n = 0, nL = rhs->length(); n < nL; ++n)
    {
      Selector_Obj r = (*rhs)[n];
      if (Wrapped_Selector_Obj wrapped = Cast<Wrapped_Selector>(r)) {
        if (wrapped->name() == ":not") {
          if (Selector_List_Obj ls = Cast<Selector_List>(wrapped->selector())) {
            ls->remove_parent_selectors(); // unverified
            if (is_superselector_of(ls, wrapped->name())) return false;
          }
        }
        if (wrapped->name() == ":matches" || wrapped->name() == ":-moz-any") {
          if (!wrapping.empty()) {
            if (wrapping != wrapped->name()) return false;
          }
          if (Selector_List_Obj ls = Cast<Selector_List>(wrapped->selector())) {
            ls->remove_parent_selectors(); // unverified
            return (is_superselector_of(ls, wrapped->name()));
          }
        }
      }
      rset.insert(r);
    }

    return contains_all(rset, lset);
  }

  bool Compound_Selector::is_universal() const
  {
    return length() == 1 && (*this)[0]->is_universal();
  }

  // create complex selector (ancestor of) from compound selector
  Complex_Selector_Obj Compound_Selector::to_complex()
  {
    // create an intermediate complex selector
    return SASS_MEMORY_NEW(Complex_Selector,
                           pstate(),
                           Complex_Selector::ANCESTOR_OF,
                           this,
                           {});
  }

  Simple_Selector_Ptr Compound_Selector::base() const {
    if (length() == 0) return 0;
    // ToDo: why is this needed?
    if (Cast<Type_Selector>((*this)[0]))
      return (*this)[0];
    return 0;
  }

  size_t Compound_Selector::hash() const
  {
    if (Selector::hash_ == 0) {
      hash_combine(Selector::hash_, std::hash<int>()(SELECTOR));
      if (length()) hash_combine(Selector::hash_, Vectorized::hash());
    }
    return Selector::hash_;
  }

  unsigned long Compound_Selector::specificity() const
  {
    int sum = 0;
    for (size_t i = 0, L = length(); i < L; ++i)
    { sum += (*this)[i]->specificity(); }
    return sum;
  }

  bool Compound_Selector::has_placeholder()
  {
    if (length() == 0) return false;
    if (Simple_Selector_Obj ss = elements().front()) {
      if (ss->has_placeholder()) return true;
    }
    return false;
  }

  bool Compound_Selector::is_empty_reference()
  {
    return length() == 1 &&
            Cast<Parent_Selector>((*this)[0]);
  }

  void Compound_Selector::append(Simple_Selector_Obj element)
  {
    Vectorized<Simple_Selector_Obj>::append(element);
    pstate_.offset += element->pstate().offset;
  }

  Compound_Selector_Ptr Compound_Selector::minus(Compound_Selector_Ptr rhs)
  {
    Compound_Selector_Ptr result = SASS_MEMORY_NEW(Compound_Selector, pstate());
    // result->has_parent_reference(has_parent_reference());

    // not very efficient because it needs to preserve order
    for (size_t i = 0, L = length(); i < L; ++i)
    {
      bool found = false;
      for (size_t j = 0, M = rhs->length(); j < M; ++j)
      {
        if (*get(i) == *rhs->get(j))
        {
          found = true;
          break;
        }
      }
      if (!found) result->append(get(i));
    }

    return result;
  }

  void Compound_Selector::mergeSources(ComplexSelectorSet& sources)
  {
    for (ComplexSelectorSet::iterator iterator = sources.begin(), endIterator = sources.end(); iterator != endIterator; ++iterator) {
      this->sources_.insert(SASS_MEMORY_CLONE(*iterator));
    }
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Complex_Selector::Complex_Selector(ParserState pstate,
                    Combinator c,
                    Compound_Selector_Obj h,
                    Complex_Selector_Obj t,
                    String_Obj r)
  : Selector(pstate),
    combinator_(c),
    head_(h), tail_(t),
    reference_(r)
  {}
  Complex_Selector::Complex_Selector(const Complex_Selector* ptr)
  : Selector(ptr),
    combinator_(ptr->combinator_),
    head_(ptr->head_), tail_(ptr->tail_),
    reference_(ptr->reference_)
  {}

  bool Complex_Selector::empty() const {
    return (!tail() || tail()->empty())
      && (!head() || head()->empty())
      && combinator_ == ANCESTOR_OF;
  }

  Complex_Selector_Obj Complex_Selector::skip_empty_reference()
  {
    if ((!head_ || !head_->length() || head_->is_empty_reference()) &&
        combinator() == Combinator::ANCESTOR_OF)
    {
      if (!tail_) return {};
      tail_->has_line_feed_ = this->has_line_feed_;
      // tail_->has_line_break_ = this->has_line_break_;
      return tail_->skip_empty_reference();
    }
    return this;
  }

  bool Complex_Selector::is_empty_ancestor() const
  {
    return (!head() || head()->length() == 0) &&
            combinator() == Combinator::ANCESTOR_OF;
  }

  size_t Complex_Selector::hash() const
  {
    if (hash_ == 0) {
      if (head_) {
        hash_combine(hash_, head_->hash());
      } else {
        hash_combine(hash_, std::hash<int>()(SELECTOR));
      }
      if (tail_) hash_combine(hash_, tail_->hash());
      if (combinator_ != ANCESTOR_OF) hash_combine(hash_, std::hash<int>()(combinator_));
    }
    return hash_;
  }

  unsigned long Complex_Selector::specificity() const
  {
    int sum = 0;
    if (head()) sum += head()->specificity();
    if (tail()) sum += tail()->specificity();
    return sum;
  }

  void Complex_Selector::set_media_block(Media_Block_Ptr mb) {
    media_block(mb);
    if (tail_) tail_->set_media_block(mb);
    if (head_) head_->set_media_block(mb);
  }

  bool Complex_Selector::has_placeholder() {
    if (head_ && head_->has_placeholder()) return true;
    if (tail_ && tail_->has_placeholder()) return true;
    return false;
  }

  const ComplexSelectorSet Complex_Selector::sources()
  {
    //s = Set.new
    //seq.map {|sseq_or_op| s.merge sseq_or_op.sources if sseq_or_op.is_a?(SimpleSequence)}
    //s

    ComplexSelectorSet srcs;

    Compound_Selector_Obj pHead = head();
    Complex_Selector_Obj  pTail = tail();

    if (pHead) {
      const ComplexSelectorSet& headSources = pHead->sources();
      srcs.insert(headSources.begin(), headSources.end());
    }

    if (pTail) {
      const ComplexSelectorSet& tailSources = pTail->sources();
      srcs.insert(tailSources.begin(), tailSources.end());
    }

    return srcs;
  }

  void Complex_Selector::addSources(ComplexSelectorSet& sources)
  {
    // members.map! {|m| m.is_a?(SimpleSequence) ? m.with_more_sources(sources) : m}
    Complex_Selector_Ptr pIter = this;
    while (pIter) {
      Compound_Selector_Ptr pHead = pIter->head();

      if (pHead) {
        pHead->mergeSources(sources);
      }

      pIter = pIter->tail();
    }
  }

  void Complex_Selector::clearSources()
  {
    Complex_Selector_Ptr pIter = this;
    while (pIter) {
      Compound_Selector_Ptr pHead = pIter->head();

      if (pHead) {
        pHead->clearSources();
      }

      pIter = pIter->tail();
    }
  }

  bool Complex_Selector::find ( bool (*f)(AST_Node_Obj) )
  {
    // check children first
    if (head_ && head_->find(f)) return true;
    if (tail_ && tail_->find(f)) return true;
    // execute last
    return f(this);
  }

  bool Complex_Selector::has_parent_ref() const
  {
    return (head() && head()->has_parent_ref()) ||
           (tail() && tail()->has_parent_ref());
  }

  bool Complex_Selector::has_real_parent_ref() const
  {
    return (head() && head()->has_real_parent_ref()) ||
           (tail() && tail()->has_real_parent_ref());
  }

  bool Complex_Selector::is_superselector_of(Compound_Selector_Ptr_Const rhs, std::string wrapping) const
  {
    return last()->head() && last()->head()->is_superselector_of(rhs, wrapping);
  }

  bool Complex_Selector::is_superselector_of(Complex_Selector_Ptr_Const rhs, std::string wrapping) const
  {
    Complex_Selector_Ptr_Const lhs = this;
    // check for selectors with leading or trailing combinators
    if (!lhs->head() || !rhs->head())
    { return false; }
    Complex_Selector_Ptr_Const l_innermost = lhs->last();
    if (l_innermost->combinator() != Complex_Selector::ANCESTOR_OF)
    { return false; }
    Complex_Selector_Ptr_Const r_innermost = rhs->last();
    if (r_innermost->combinator() != Complex_Selector::ANCESTOR_OF)
    { return false; }
    // more complex (i.e., longer) selectors are always more specific
    size_t l_len = lhs->length(), r_len = rhs->length();
    if (l_len > r_len)
    { return false; }

    if (l_len == 1)
    { return lhs->head()->is_superselector_of(rhs->last()->head(), wrapping); }

    // we have to look one tail deeper, since we cary the
    // combinator around for it (which is important here)
    if (rhs->tail() && lhs->tail() && combinator() != Complex_Selector::ANCESTOR_OF) {
      Complex_Selector_Obj lhs_tail = lhs->tail();
      Complex_Selector_Obj rhs_tail = rhs->tail();
      if (lhs_tail->combinator() != rhs_tail->combinator()) return false;
      if (lhs_tail->head() && !rhs_tail->head()) return false;
      if (!lhs_tail->head() && rhs_tail->head()) return false;
      if (lhs_tail->head() && rhs_tail->head()) {
        if (!lhs_tail->head()->is_superselector_of(rhs_tail->head())) return false;
      }
    }

    bool found = false;
    Complex_Selector_Ptr_Const marker = rhs;
    for (size_t i = 0, L = rhs->length(); i < L; ++i) {
      if (i == L-1)
      { return false; }
      if (lhs->head() && marker->head() && lhs->head()->is_superselector_of(marker->head(), wrapping))
      { found = true; break; }
      marker = marker->tail();
    }
    if (!found)
    { return false; }

    /*
      Hmm, I hope I have the logic right:

      if lhs has a combinator:
        if !(marker has a combinator) return false
        if !(lhs.combinator == '~' ? marker.combinator != '>' : lhs.combinator == marker.combinator) return false
        return lhs.tail-without-innermost.is_superselector_of(marker.tail-without-innermost)
      else if marker has a combinator:
        if !(marker.combinator == ">") return false
        return lhs.tail.is_superselector_of(marker.tail)
      else
        return lhs.tail.is_superselector_of(marker.tail)
    */
    if (lhs->combinator() != Complex_Selector::ANCESTOR_OF)
    {
      if (marker->combinator() == Complex_Selector::ANCESTOR_OF)
      { return false; }
      if (!(lhs->combinator() == Complex_Selector::PRECEDES ? marker->combinator() != Complex_Selector::PARENT_OF : lhs->combinator() == marker->combinator()))
      { return false; }
      return lhs->tail()->is_superselector_of(marker->tail());
    }
    else if (marker->combinator() != Complex_Selector::ANCESTOR_OF)
    {
      if (marker->combinator() != Complex_Selector::PARENT_OF)
      { return false; }
      return lhs->tail()->is_superselector_of(marker->tail());
    }
    return lhs->tail()->is_superselector_of(marker->tail());
  }

  size_t Complex_Selector::length() const
  {
    // TODO: make this iterative
    if (!tail()) return 1;
    return 1 + tail()->length();
  }

  // append another complex selector at the end
  // check if we need to append some headers
  // then we need to check for the combinator
  // only then we can safely set the new tail
  void Complex_Selector::append(Complex_Selector_Obj ss, Backtraces& traces)
  {

    Complex_Selector_Obj t = ss->tail();
    Combinator c = ss->combinator();
    String_Obj r = ss->reference();
    Compound_Selector_Obj h = ss->head();

    if (ss->has_line_feed()) has_line_feed(true);
    if (ss->has_line_break()) has_line_break(true);

    // append old headers
    if (h && h->length()) {
      if (last()->combinator() != ANCESTOR_OF && c != ANCESTOR_OF) {
        traces.push_back(Backtrace(pstate()));
        throw Exception::InvalidParent(this, traces, ss);
      } else if (last()->head_ && last()->head_->length()) {
        Compound_Selector_Obj rh = last()->head();
        size_t i;
        size_t L = h->length();
        if (Cast<Type_Selector>(h->first())) {
          if (Class_Selector_Ptr cs = Cast<Class_Selector>(rh->last())) {
            Class_Selector_Ptr sqs = SASS_MEMORY_COPY(cs);
            sqs->name(sqs->name() + (*h)[0]->name());
            sqs->pstate((*h)[0]->pstate());
            (*rh)[rh->length()-1] = sqs;
            rh->pstate(h->pstate());
            for (i = 1; i < L; ++i) rh->append((*h)[i]);
          } else if (Id_Selector_Ptr is = Cast<Id_Selector>(rh->last())) {
            Id_Selector_Ptr sqs = SASS_MEMORY_COPY(is);
            sqs->name(sqs->name() + (*h)[0]->name());
            sqs->pstate((*h)[0]->pstate());
            (*rh)[rh->length()-1] = sqs;
            rh->pstate(h->pstate());
            for (i = 1; i < L; ++i) rh->append((*h)[i]);
          } else if (Type_Selector_Ptr ts = Cast<Type_Selector>(rh->last())) {
            Type_Selector_Ptr tss = SASS_MEMORY_COPY(ts);
            tss->name(tss->name() + (*h)[0]->name());
            tss->pstate((*h)[0]->pstate());
            (*rh)[rh->length()-1] = tss;
            rh->pstate(h->pstate());
            for (i = 1; i < L; ++i) rh->append((*h)[i]);
          } else if (Placeholder_Selector_Ptr ps = Cast<Placeholder_Selector>(rh->last())) {
            Placeholder_Selector_Ptr pss = SASS_MEMORY_COPY(ps);
            pss->name(pss->name() + (*h)[0]->name());
            pss->pstate((*h)[0]->pstate());
            (*rh)[rh->length()-1] = pss;
            rh->pstate(h->pstate());
            for (i = 1; i < L; ++i) rh->append((*h)[i]);
          } else {
            last()->head_->concat(h);
          }
        } else {
          last()->head_->concat(h);
        }
      } else if (last()->head_) {
        last()->head_->concat(h);
      }
    } else {
      // std::cerr << "has no or empty head\n";
    }

    Complex_Selector_Ptr last = mutable_last();
    if (last) {
      if (last->combinator() != ANCESTOR_OF && c != ANCESTOR_OF) {
        Complex_Selector_Ptr inter = SASS_MEMORY_NEW(Complex_Selector, pstate());
        inter->reference(r);
        inter->combinator(c);
        inter->tail(t);
        last->tail(inter);
      } else {
        if (last->combinator() == ANCESTOR_OF) {
          last->combinator(c);
          last->reference(r);
        }
        last->tail(t);
      }
    }

  }

  Selector_List_Ptr Complex_Selector::resolve_parent_refs(SelectorStack& pstack, Backtraces& traces, bool implicit_parent)
  {
    Complex_Selector_Obj tail = this->tail();
    Compound_Selector_Obj head = this->head();
    Selector_List_Ptr parents = pstack.back();

    if (!this->has_real_parent_ref() && !implicit_parent) {
      Selector_List_Ptr retval = SASS_MEMORY_NEW(Selector_List, pstate(), 1);
      retval->append(this);
      return retval;
    }

    // first resolve_parent_refs the tail (which may return an expanded list)
    Selector_List_Obj tails = tail ? tail->resolve_parent_refs(pstack, traces, implicit_parent) : 0;

    if (head && head->length() > 0) {

      Selector_List_Obj retval;
      // we have a parent selector in a simple compound list
      // mix parent complex selector into the compound list
      if (Cast<Parent_Selector>((*head)[0])) {
        retval = SASS_MEMORY_NEW(Selector_List, pstate());

        // it turns out that real parent references reach
        // across @at-root rules, which comes unexpected
        if (parents == NULL && head->has_real_parent_ref()) {
          int i = pstack.size() - 1;
          while (!parents && i > -1) {
            parents = pstack.at(i--);
          }
        }

        if (parents && parents->length()) {
          if (tails && tails->length() > 0) {
            for (size_t n = 0, nL = tails->length(); n < nL; ++n) {
              for (size_t i = 0, iL = parents->length(); i < iL; ++i) {
                Complex_Selector_Obj t = (*tails)[n];
                Complex_Selector_Obj parent = (*parents)[i];
                Complex_Selector_Obj s = SASS_MEMORY_CLONE(parent);
                Complex_Selector_Obj ss = SASS_MEMORY_CLONE(this);
                ss->tail(t ? SASS_MEMORY_CLONE(t) : NULL);
                Compound_Selector_Obj h = SASS_MEMORY_COPY(head_);
                // remove parent selector from sequence
                if (h->length()) {
                  h->erase(h->begin());
                  ss->head(h);
                } else {
                  ss->head({});
                }
                // adjust for parent selector (1 char)
                // if (h->length()) {
                //   ParserState state(h->at(0)->pstate());
                //   state.offset.column += 1;
                //   state.column -= 1;
                //   (*h)[0]->pstate(state);
                // }
                // keep old parser state
                s->pstate(pstate());
                // append new tail
                s->append(ss, traces);
                retval->append(s);
              }
            }
          }
          // have no tails but parents
          // loop above is inside out
          else {
            for (size_t i = 0, iL = parents->length(); i < iL; ++i) {
              Complex_Selector_Obj parent = (*parents)[i];
              Complex_Selector_Obj s = SASS_MEMORY_CLONE(parent);
              Complex_Selector_Obj ss = SASS_MEMORY_CLONE(this);
              // this is only if valid if the parent has no trailing op
              // otherwise we cannot append more simple selectors to head
              if (parent->last()->combinator() != ANCESTOR_OF) {
                traces.push_back(Backtrace(pstate()));
                throw Exception::InvalidParent(parent, traces, ss);
              }
              ss->tail(tail ? SASS_MEMORY_CLONE(tail) : NULL);
              Compound_Selector_Obj h = SASS_MEMORY_COPY(head_);
              // remove parent selector from sequence
              if (h->length()) {
                h->erase(h->begin());
                ss->head(h);
              } else {
                ss->head({});
              }
              // \/ IMO ruby sass bug \/
              ss->has_line_feed(false);
              // adjust for parent selector (1 char)
              // if (h->length()) {
              //   ParserState state(h->at(0)->pstate());
              //   state.offset.column += 1;
              //   state.column -= 1;
              //   (*h)[0]->pstate(state);
              // }
              // keep old parser state
              s->pstate(pstate());
              // append new tail
              s->append(ss, traces);
              retval->append(s);
            }
          }
        }
        // have no parent but some tails
        else {
          if (tails && tails->length() > 0) {
            for (size_t n = 0, nL = tails->length(); n < nL; ++n) {
              Complex_Selector_Obj cpy = SASS_MEMORY_CLONE(this);
              cpy->tail(SASS_MEMORY_CLONE(tails->at(n)));
              cpy->head(SASS_MEMORY_NEW(Compound_Selector, head->pstate()));
              for (size_t i = 1, L = this->head()->length(); i < L; ++i)
                cpy->head()->append((*this->head())[i]);
              if (!cpy->head()->length()) cpy->head({});
              retval->append(cpy->skip_empty_reference());
            }
          }
          // have no parent nor tails
          else {
            Complex_Selector_Obj cpy = SASS_MEMORY_CLONE(this);
            cpy->head(SASS_MEMORY_NEW(Compound_Selector, head->pstate()));
            for (size_t i = 1, L = this->head()->length(); i < L; ++i)
              cpy->head()->append((*this->head())[i]);
            if (!cpy->head()->length()) cpy->head({});
            retval->append(cpy->skip_empty_reference());
          }
        }
      }
      // no parent selector in head
      else {
        retval = this->tails(tails);
      }

      for (Simple_Selector_Obj ss : head->elements()) {
        if (Wrapped_Selector_Ptr ws = Cast<Wrapped_Selector>(ss)) {
          if (Selector_List_Ptr sl = Cast<Selector_List>(ws->selector())) {
            if (parents) ws->selector(sl->resolve_parent_refs(pstack, traces, implicit_parent));
          }
        }
      }

      return retval.detach();

    }
    // has no head
    return this->tails(tails);
  }

  Selector_List_Ptr Complex_Selector::tails(Selector_List_Ptr tails)
  {
    Selector_List_Ptr rv = SASS_MEMORY_NEW(Selector_List, pstate_);
    if (tails && tails->length()) {
      for (size_t i = 0, iL = tails->length(); i < iL; ++i) {
        Complex_Selector_Obj pr = SASS_MEMORY_CLONE(this);
        pr->tail(tails->at(i));
        rv->append(pr);
      }
    }
    else {
      rv->append(this);
    }
    return rv;
  }

  // return the last tail that is defined
  Complex_Selector_Ptr_Const Complex_Selector::first() const
  {
    // declare variables used in loop
    Complex_Selector_Ptr_Const cur = this;
    Compound_Selector_Ptr_Const head;
    // processing loop
    while (cur)
    {
      // get the head
      head = cur->head_.ptr();
      // abort (and return) if it is not a parent selector
      if (!head || head->length() != 1 || !Cast<Parent_Selector>((*head)[0])) {
        break;
      }
      // advance to next
      cur = cur->tail_;
    }
    // result
    return cur;
  }

  Complex_Selector_Ptr Complex_Selector::mutable_first()
  {
    return const_cast<Complex_Selector_Ptr>(first());
  }

  // return the last tail that is defined
  Complex_Selector_Ptr_Const Complex_Selector::last() const
  {
    Complex_Selector_Ptr_Const cur = this;
    Complex_Selector_Ptr_Const nxt = cur;
    // loop until last
    while (nxt) {
      cur = nxt;
      nxt = cur->tail_.ptr();
    }
    return cur;
  }

  Complex_Selector_Ptr Complex_Selector::mutable_last()
  {
    return const_cast<Complex_Selector_Ptr>(last());
  }

  Complex_Selector::Combinator Complex_Selector::clear_innermost()
  {
    Combinator c;
    if (!tail() || tail()->tail() == nullptr)
    { c = combinator(); combinator(ANCESTOR_OF); tail({}); }
    else
    { c = tail_->clear_innermost(); }
    return c;
  }

  void Complex_Selector::set_innermost(Complex_Selector_Obj val, Combinator c)
  {
    if (!tail_)
    { tail_ = val; combinator(c); }
    else
    { tail_->set_innermost(val, c); }
  }

  void Complex_Selector::cloneChildren()
  {
    if (head()) head(SASS_MEMORY_CLONE(head()));
    if (tail()) tail(SASS_MEMORY_CLONE(tail()));
  }

  // it's a superselector if every selector of the right side
  // list is a superselector of the given left side selector
  bool Complex_Selector::is_superselector_of(Selector_List_Ptr_Const sub, std::string wrapping) const
  {
    // Check every rhs selector against left hand list
    for(size_t i = 0, L = sub->length(); i < L; ++i) {
      if (!is_superselector_of((*sub)[i], wrapping)) return false;
    }
    return true;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  Selector_List::Selector_List(ParserState pstate, size_t s)
  : Selector(pstate),
    Vectorized<Complex_Selector_Obj>(s),
    schema_({}),
    wspace_(0)
  { }
  Selector_List::Selector_List(const Selector_List* ptr)
  : Selector(ptr),
    Vectorized<Complex_Selector_Obj>(*ptr),
    schema_(ptr->schema_),
    wspace_(ptr->wspace_)
  { }

  bool Selector_List::find ( bool (*f)(AST_Node_Obj) )
  {
    // check children first
    for (Complex_Selector_Obj sel : elements()) {
      if (sel->find(f)) return true;
    }
    // execute last
    return f(this);
  }

  Selector_List_Obj Selector_List::eval(Eval& eval)
  {
    Selector_List_Obj list = schema() ?
      eval(schema()) : eval(this);
    list->schema(schema());
    return list;
  }

  Selector_List_Ptr Selector_List::resolve_parent_refs(SelectorStack& pstack, Backtraces& traces, bool implicit_parent)
  {
    if (!this->has_parent_ref()) return this;
    Selector_List_Ptr ss = SASS_MEMORY_NEW(Selector_List, pstate());
    for (size_t si = 0, sL = this->length(); si < sL; ++si) {
      Selector_List_Obj rv = at(si)->resolve_parent_refs(pstack, traces, implicit_parent);
      ss->concat(rv);
    }
    return ss;
  }

  void Selector_List::cloneChildren()
  {
    for (size_t i = 0, l = length(); i < l; i++) {
      at(i) = SASS_MEMORY_CLONE(at(i));
    }
  }

  // remove parent selector references
  // basically unwraps parsed selectors
  void Selector_List::remove_parent_selectors()
  {
    // Check every rhs selector against left hand list
    for(size_t i = 0, L = length(); i < L; ++i) {
      if (!(*this)[i]->head()) continue;
      if ((*this)[i]->head()->is_empty_reference()) {
        // simply move to the next tail if we have "no" combinator
        if ((*this)[i]->combinator() == Complex_Selector::ANCESTOR_OF) {
          if ((*this)[i]->tail()) {
            if ((*this)[i]->has_line_feed()) {
              (*this)[i]->tail()->has_line_feed(true);
            }
            (*this)[i] = (*this)[i]->tail();
          }
        }
        // otherwise remove the first item from head
        else {
          (*this)[i]->head()->erase((*this)[i]->head()->begin());
        }
      }
    }
  }

  bool Selector_List::has_parent_ref() const
  {
    for (Complex_Selector_Obj s : elements()) {
      if (s && s->has_parent_ref()) return true;
    }
    return false;
  }

  bool Selector_List::has_real_parent_ref() const
  {
    for (Complex_Selector_Obj s : elements()) {
      if (s && s->has_real_parent_ref()) return true;
    }
    return false;
  }

  void Selector_List::adjust_after_pushing(Complex_Selector_Obj c)
  {
    // if (c->has_reference())   has_reference(true);
  }

  // it's a superselector if every selector of the right side
  // list is a superselector of the given left side selector
  bool Selector_List::is_superselector_of(Selector_List_Ptr_Const sub, std::string wrapping) const
  {
    // Check every rhs selector against left hand list
    for(size_t i = 0, L = sub->length(); i < L; ++i) {
      if (!is_superselector_of((*sub)[i], wrapping)) return false;
    }
    return true;
  }

  // it's a superselector if every selector on the right side
  // is a superselector of any one of the left side selectors
  bool Selector_List::is_superselector_of(Compound_Selector_Ptr_Const sub, std::string wrapping) const
  {
    // Check every lhs selector against right hand
    for(size_t i = 0, L = length(); i < L; ++i) {
      if ((*this)[i]->is_superselector_of(sub, wrapping)) return true;
    }
    return false;
  }

  // it's a superselector if every selector on the right side
  // is a superselector of any one of the left side selectors
  bool Selector_List::is_superselector_of(Complex_Selector_Ptr_Const sub, std::string wrapping) const
  {
    // Check every lhs selector against right hand
    for(size_t i = 0, L = length(); i < L; ++i) {
      if ((*this)[i]->is_superselector_of(sub)) return true;
    }
    return false;
  }

  void Selector_List::populate_extends(Selector_List_Obj extendee, Subset_Map& extends)
  {

    Selector_List_Ptr extender = this;
    for (auto complex_sel : extendee->elements()) {
      Complex_Selector_Obj c = complex_sel;


      // Ignore any parent selectors, until we find the first non Selectorerence head
      Compound_Selector_Obj compound_sel = c->head();
      Complex_Selector_Obj pIter = complex_sel;
      while (pIter) {
        Compound_Selector_Obj pHead = pIter->head();
        if (pHead && Cast<Parent_Selector>(pHead->elements()[0]) == NULL) {
          compound_sel = pHead;
          break;
        }

        pIter = pIter->tail();
      }

      if (!pIter->head() || pIter->tail()) {
        coreError("nested selectors may not be extended", c->pstate());
      }

      compound_sel->is_optional(extendee->is_optional());

      for (size_t i = 0, L = extender->length(); i < L; ++i) {
        extends.put(compound_sel, std::make_pair((*extender)[i], compound_sel));
      }
    }
  };

  size_t Selector_List::hash() const
  {
    if (Selector::hash_ == 0) {
      if (empty()) {
        hash_combine(Selector::hash_, std::hash<int>()(SELECTOR));
      } else {
        hash_combine(Selector::hash_, Vectorized::hash());
      }
    }
    return Selector::hash_;
  }

  unsigned long Selector_List::specificity() const
  {
    unsigned long sum = 0;
    unsigned long specificity;
    for (size_t i = 0, L = length(); i < L; ++i)
    {
      specificity = (*this)[i]->specificity();
      if (sum < specificity) sum = specificity;
    }
    return sum;
  }

  void Selector_List::set_media_block(Media_Block_Ptr mb)
  {
    media_block(mb);
    for (Complex_Selector_Obj cs : elements()) {
      cs->set_media_block(mb);
    }
  }

  bool Selector_List::has_placeholder()
  {
    for (Complex_Selector_Obj cs : elements()) {
      if (cs->has_placeholder()) return true;
    }
    return false;
  }

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////

  IMPLEMENT_AST_OPERATORS(Selector_Schema);
  IMPLEMENT_AST_OPERATORS(Placeholder_Selector);
  IMPLEMENT_AST_OPERATORS(Parent_Selector);
  IMPLEMENT_AST_OPERATORS(Attribute_Selector);
  IMPLEMENT_AST_OPERATORS(Compound_Selector);
  IMPLEMENT_AST_OPERATORS(Complex_Selector);
  IMPLEMENT_AST_OPERATORS(Type_Selector);
  IMPLEMENT_AST_OPERATORS(Class_Selector);
  IMPLEMENT_AST_OPERATORS(Id_Selector);
  IMPLEMENT_AST_OPERATORS(Pseudo_Selector);
  IMPLEMENT_AST_OPERATORS(Wrapped_Selector);
  IMPLEMENT_AST_OPERATORS(Selector_List);

}

namespace Sass {

  namespace Functions {

    /////////////////
    // MAP FUNCTIONS
    /////////////////

    Signature map_get_sig = "map-get($map, $key)";
    BUILT_IN(map_get)
    {
      // leaks for "map-get((), foo)" if not Obj
      // investigate why this is (unexpected)
      Map_Obj m = ARGM("$map", Map);
      Expression_Obj v = ARG("$key", Expression);
      try {
        Value_Obj val = m->at(v);
        if (!val) return SASS_MEMORY_NEW(Null, pstate);
        val->set_delayed(false);
        return val.detach();
      } catch (const std::out_of_range&) {
        return SASS_MEMORY_NEW(Null, pstate);
      }
      catch (...) { throw; }
    }

    Signature map_has_key_sig = "map-has-key($map, $key)";
    BUILT_IN(map_has_key)
    {
      Map_Obj m = ARGM("$map", Map);
      Expression_Obj v = ARG("$key", Expression);
      return SASS_MEMORY_NEW(Boolean, pstate, m->has(v));
    }

    Signature map_keys_sig = "map-keys($map)";
    BUILT_IN(map_keys)
    {
      Map_Obj m = ARGM("$map", Map);
      List_Ptr result = SASS_MEMORY_NEW(List, pstate, m->length(), SASS_COMMA);
      for ( auto key : m->keys()) {
        result->append(key);
      }
      return result;
    }

    Signature map_values_sig = "map-values($map)";
    BUILT_IN(map_values)
    {
      Map_Obj m = ARGM("$map", Map);
      List_Ptr result = SASS_MEMORY_NEW(List, pstate, m->length(), SASS_COMMA);
      for ( auto key : m->keys()) {
        result->append(m->at(key));
      }
      return result;
    }

    Signature map_merge_sig = "map-merge($map1, $map2)";
    BUILT_IN(map_merge)
    {
      Map_Obj m1 = ARGM("$map1", Map);
      Map_Obj m2 = ARGM("$map2", Map);

      size_t len = m1->length() + m2->length();
      Map_Ptr result = SASS_MEMORY_NEW(Map, pstate, len);
      // concat not implemented for maps
      *result += m1;
      *result += m2;
      return result;
    }

    Signature map_remove_sig = "map-remove($map, $keys...)";
    BUILT_IN(map_remove)
    {
      bool remove;
      Map_Obj m = ARGM("$map", Map);
      List_Obj arglist = ARG("$keys", List);
      Map_Ptr result = SASS_MEMORY_NEW(Map, pstate, 1);
      for (auto key : m->keys()) {
        remove = false;
        for (size_t j = 0, K = arglist->length(); j < K && !remove; ++j) {
          remove = Operators::eq(key, arglist->value_at_index(j));
        }
        if (!remove) *result << std::make_pair(key, m->at(key));
      }
      return result;
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {


  Offset::Offset(const char chr)
  : line(chr == '\n' ? 1 : 0),
    column(chr == '\n' ? 0 : 1)
  {}

  Offset::Offset(const char* string)
  : line(0), column(0)
  {
    *this = inc(string, string + strlen(string));
  }

  Offset::Offset(const std::string& text)
  : line(0), column(0)
  {
    *this = inc(text.c_str(), text.c_str() + text.size());
  }

  Offset::Offset(const size_t line, const size_t column)
  : line(line), column(column) { }

  // init/create instance from const char substring
  Offset Offset::init(const char* beg, const char* end)
  {
    Offset offset(0, 0);
    if (end == 0) {
      end += strlen(beg);
    }
    offset.add(beg, end);
    return offset;
  }

  // increase offset by given string (mostly called by lexer)
  // increase line counter and count columns on the last line
  Offset Offset::add(const char* begin, const char* end)
  {
    if (end == 0) return *this;
    while (begin < end && *begin) {
      if (*begin == '\n') {
        ++ line;
        // start new line
        column = 0;
      } else {
        // do not count any utf8 continuation bytes
        // https://stackoverflow.com/a/9356203/1550314
        // https://en.wikipedia.org/wiki/UTF-8#Description
        unsigned char chr = *begin;
        // skip over 10xxxxxx
        // is 1st bit not set
        if ((chr & 128) == 0) {
          // regular ascii char
          column += 1;
        }
        // is 2nd bit not set
        else if ((chr & 64) == 0) {
          // first utf8 byte
          column += 1;
        }
      }
      ++ begin;
    }
    return *this;
  }

  // increase offset by given string (mostly called by lexer)
  // increase line counter and count columns on the last line
  Offset Offset::inc(const char* begin, const char* end) const
  {
    Offset offset(line, column);
    offset.add(begin, end);
    return offset;
  }

  bool Offset::operator== (const Offset &pos) const
  {
    return line == pos.line && column == pos.column;
  }

  bool Offset::operator!= (const Offset &pos) const
  {
    return line != pos.line || column != pos.column;
  }

  void Offset::operator+= (const Offset &off)
  {
    *this = Offset(line + off.line, off.line > 0 ? off.column : column + off.column);
  }

  Offset Offset::operator+ (const Offset &off) const
  {
    return Offset(line + off.line, off.line > 0 ? off.column : column + off.column);
  }

  Offset Offset::operator- (const Offset &off) const
  {
    return Offset(line - off.line, off.line == line ? column - off.column : column);
  }

  Position::Position(const size_t file)
  : Offset(0, 0), file(file) { }

  Position::Position(const size_t file, const Offset& offset)
  : Offset(offset), file(file) { }

  Position::Position(const size_t line, const size_t column)
  : Offset(line, column), file(-1) { }

  Position::Position(const size_t file, const size_t line, const size_t column)
  : Offset(line, column), file(file) { }


  ParserState::ParserState(const char* path, const char* src, const size_t file)
  : Position(file, 0, 0), path(path), src(src), offset(0, 0), token() { }

  ParserState::ParserState(const char* path, const char* src, const Position& position, Offset offset)
  : Position(position), path(path), src(src), offset(offset), token() { }

  ParserState::ParserState(const char* path, const char* src, const Token& token, const Position& position, Offset offset)
  : Position(position), path(path), src(src), offset(offset), token(token) { }

  Position Position::add(const char* begin, const char* end)
  {
    Offset::add(begin, end);
    return *this;
  }

  Position Position::inc(const char* begin, const char* end) const
  {
    Offset offset(line, column);
    offset = offset.inc(begin, end);
    return Position(file, offset);
  }

  bool Position::operator== (const Position &pos) const
  {
    return file == pos.file && line == pos.line && column == pos.column;
  }

  bool Position::operator!= (const Position &pos) const
  {
    return file == pos.file || line != pos.line || column != pos.column;
  }

  void Position::operator+= (const Offset &off)
  {
    *this = Position(file, line + off.line, off.line > 0 ? off.column : column + off.column);
  }

  const Position Position::operator+ (const Offset &off) const
  {
    return Position(file, line + off.line, off.line > 0 ? off.column : column + off.column);
  }

  const Offset Position::operator- (const Offset &off) const
  {
    return Offset(line - off.line, off.line == line ? column - off.column : column);
  }

  /* not used anymore - remove?
  std::ostream& operator<<(std::ostream& strm, const Offset& off)
  {
    if (off.line == string::npos) strm << "-1:"; else strm << off.line << ":";
    if (off.column == string::npos) strm << "-1"; else strm << off.column;
    return strm;
  } */

  /* not used anymore - remove?
  std::ostream& operator<<(std::ostream& strm, const Position& pos)
  {
    if (pos.file != string::npos) strm << pos.file << ":";
    if (pos.line == string::npos) strm << "-1:"; else strm << pos.line << ":";
    if (pos.column == string::npos) strm << "-1"; else strm << pos.column;
    return strm;
  } */

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


namespace Sass {

  union Sass_Value* AST2C::operator()(Boolean_Ptr b)
  { return sass_make_boolean(b->value()); }

  union Sass_Value* AST2C::operator()(Number_Ptr n)
  { return sass_make_number(n->value(), n->unit().c_str()); }

  union Sass_Value* AST2C::operator()(Custom_Warning_Ptr w)
  { return sass_make_warning(w->message().c_str()); }

  union Sass_Value* AST2C::operator()(Custom_Error_Ptr e)
  { return sass_make_error(e->message().c_str()); }

  union Sass_Value* AST2C::operator()(Color_RGBA_Ptr c)
  { return sass_make_color(c->r(), c->g(), c->b(), c->a()); }

  union Sass_Value* AST2C::operator()(Color_HSLA_Ptr c)
  { return operator()(c->toRGBA()); }

  union Sass_Value* AST2C::operator()(String_Constant_Ptr s)
  {
    if (s->quote_mark()) {
      return sass_make_qstring(s->value().c_str());
    } else {
      return sass_make_string(s->value().c_str());
    }
  }

  union Sass_Value* AST2C::operator()(String_Quoted_Ptr s)
  { return sass_make_qstring(s->value().c_str()); }

  union Sass_Value* AST2C::operator()(List_Ptr l)
  {
    union Sass_Value* v = sass_make_list(l->length(), l->separator(), l->is_bracketed());
    for (size_t i = 0, L = l->length(); i < L; ++i) {
      sass_list_set_value(v, i, (*l)[i]->perform(this));
    }
    return v;
  }

  union Sass_Value* AST2C::operator()(Map_Ptr m)
  {
    union Sass_Value* v = sass_make_map(m->length());
    int i = 0;
    for (auto key : m->keys()) {
      sass_map_set_key(v, i, key->perform(this));
      sass_map_set_value(v, i, m->at(key)->perform(this));
      i++;
    }
    return v;
  }

  union Sass_Value* AST2C::operator()(Arguments_Ptr a)
  {
    union Sass_Value* v = sass_make_list(a->length(), SASS_COMMA, false);
    for (size_t i = 0, L = a->length(); i < L; ++i) {
      sass_list_set_value(v, i, (*a)[i]->perform(this));
    }
    return v;
  }

  union Sass_Value* AST2C::operator()(Argument_Ptr a)
  { return a->value()->perform(this); }

  // not strictly necessary because of the fallback
  union Sass_Value* AST2C::operator()(Null_Ptr n)
  { return sass_make_null(); }

};
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.

#include <cstdint>
#include <random>


#ifdef __MINGW32__
#include "windows.h"
#include "wincrypt.h"
#endif

namespace Sass {

  namespace Functions {

    #ifdef __MINGW32__
      uint64_t GetSeed()
      {
        HCRYPTPROV hp = 0;
        BYTE rb[8];
        CryptAcquireContext(&hp, 0, 0, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT);
        CryptGenRandom(hp, sizeof(rb), rb);
        CryptReleaseContext(hp, 0);

        uint64_t seed;
        memcpy(&seed, &rb[0], sizeof(seed));

        return seed;
      }
    #else
      uint64_t GetSeed()
      {
        std::random_device rd;
        return rd();
      }
    #endif

    // note: the performance of many implementations of
    // random_device degrades sharply once the entropy pool
    // is exhausted. For practical use, random_device is
    // generally only used to seed a PRNG such as mt19937.
    static std::mt19937 rand(static_cast<unsigned int>(GetSeed()));

    ///////////////////
    // NUMBER FUNCTIONS
    ///////////////////

    Signature percentage_sig = "percentage($number)";
    BUILT_IN(percentage)
    {
      Number_Obj n = ARGN("$number");
      if (!n->is_unitless()) error("argument $number of `" + std::string(sig) + "` must be unitless", pstate, traces);
      return SASS_MEMORY_NEW(Number, pstate, n->value() * 100, "%");
    }

    Signature round_sig = "round($number)";
    BUILT_IN(round)
    {
      Number_Obj r = ARGN("$number");
      r->value(Sass::round(r->value(), ctx.c_options.precision));
      r->pstate(pstate);
      return r.detach();
    }

    Signature ceil_sig = "ceil($number)";
    BUILT_IN(ceil)
    {
      Number_Obj r = ARGN("$number");
      r->value(std::ceil(r->value()));
      r->pstate(pstate);
      return r.detach();
    }

    Signature floor_sig = "floor($number)";
    BUILT_IN(floor)
    {
      Number_Obj r = ARGN("$number");
      r->value(std::floor(r->value()));
      r->pstate(pstate);
      return r.detach();
    }

    Signature abs_sig = "abs($number)";
    BUILT_IN(abs)
    {
      Number_Obj r = ARGN("$number");
      r->value(std::abs(r->value()));
      r->pstate(pstate);
      return r.detach();
    }

    Signature min_sig = "min($numbers...)";
    BUILT_IN(min)
    {
      List_Ptr arglist = ARG("$numbers", List);
      Number_Obj least;
      for (size_t i = 0, L = arglist->length(); i < L; ++i) {
        Expression_Obj val = arglist->value_at_index(i);
        Number_Obj xi = Cast<Number>(val);
        if (!xi) {
          error("\"" + val->to_string(ctx.c_options) + "\" is not a number for `min'", pstate, traces);
        }
        if (least) {
          if (*xi < *least) least = xi;
        } else least = xi;
      }
      return least.detach();
    }

    Signature max_sig = "max($numbers...)";
    BUILT_IN(max)
    {
      List_Ptr arglist = ARG("$numbers", List);
      Number_Obj greatest;
      for (size_t i = 0, L = arglist->length(); i < L; ++i) {
        Expression_Obj val = arglist->value_at_index(i);
        Number_Obj xi = Cast<Number>(val);
        if (!xi) {
          error("\"" + val->to_string(ctx.c_options) + "\" is not a number for `max'", pstate, traces);
        }
        if (greatest) {
          if (*greatest < *xi) greatest = xi;
        } else greatest = xi;
      }
      return greatest.detach();
    }

    Signature random_sig = "random($limit:false)";
    BUILT_IN(random)
    {
      AST_Node_Obj arg = env["$limit"];
      Value_Ptr v = Cast<Value>(arg);
      Number_Ptr l = Cast<Number>(arg);
      Boolean_Ptr b = Cast<Boolean>(arg);
      if (l) {
        double lv = l->value();
        if (lv < 1) {
          std::stringstream err;
          err << "$limit " << lv << " must be greater than or equal to 1 for `random'";
          error(err.str(), pstate, traces);
        }
        bool eq_int = std::fabs(trunc(lv) - lv) < NUMBER_EPSILON;
        if (!eq_int) {
          std::stringstream err;
          err << "Expected $limit to be an integer but got " << lv << " for `random'";
          error(err.str(), pstate, traces);
        }
        std::uniform_real_distribution<> distributor(1, lv + 1);
        uint_fast32_t distributed = static_cast<uint_fast32_t>(distributor(rand));
        return SASS_MEMORY_NEW(Number, pstate, (double)distributed);
      }
      else if (b) {
        std::uniform_real_distribution<> distributor(0, 1);
        double distributed = static_cast<double>(distributor(rand));
        return SASS_MEMORY_NEW(Number, pstate, distributed);
      } else if (v) {
        traces.push_back(Backtrace(pstate));
        throw Exception::InvalidArgumentType(pstate, traces, "random", "$limit", "number", v);
      } else {
        traces.push_back(Backtrace(pstate));
        throw Exception::InvalidArgumentType(pstate, traces, "random", "$limit", "number");
      }
    }

    Signature unique_id_sig = "unique-id()";
    BUILT_IN(unique_id)
    {
      std::stringstream ss;
      std::uniform_real_distribution<> distributor(0, 4294967296); // 16^8
      uint_fast32_t distributed = static_cast<uint_fast32_t>(distributor(rand));
      ss << "u" << std::setfill('0') << std::setw(8) << std::hex << distributed;
      return SASS_MEMORY_NEW(String_Quoted, pstate, ss.str());
    }

    Signature unit_sig = "unit($number)";
    BUILT_IN(unit)
    {
      Number_Obj arg = ARGN("$number");
      std::string str(quote(arg->unit(), '"'));
      return SASS_MEMORY_NEW(String_Quoted, pstate, str);
    }

    Signature unitless_sig = "unitless($number)";
    BUILT_IN(unitless)
    {
      Number_Obj arg = ARGN("$number");
      bool unitless = arg->is_unitless();
      return SASS_MEMORY_NEW(Boolean, pstate, unitless);
    }

    Signature comparable_sig = "comparable($number-1, $number-2)";
    BUILT_IN(comparable)
    {
      Number_Obj n1 = ARGN("$number-1");
      Number_Obj n2 = ARGN("$number-2");
      if (n1->is_unitless() || n2->is_unitless()) {
        return SASS_MEMORY_NEW(Boolean, pstate, true);
      }
      // normalize into main units
      n1->normalize(); n2->normalize();
      Units &lhs_unit = *n1, &rhs_unit = *n2;
      bool is_comparable = (lhs_unit == rhs_unit);
      return SASS_MEMORY_NEW(Boolean, pstate, is_comparable);
    }

  }

}
// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {
  using namespace Constants;

  namespace Prelexer {

    //####################################
    // BASIC CHARACTER MATCHERS
    //####################################

    // Match standard control chars
    const char* kwd_at(const char* src) { return exactly<'@'>(src); }
    const char* kwd_dot(const char* src) { return exactly<'.'>(src); }
    const char* kwd_comma(const char* src) { return exactly<','>(src); };
    const char* kwd_colon(const char* src) { return exactly<':'>(src); };
    const char* kwd_star(const char* src) { return exactly<'*'>(src); };
    const char* kwd_plus(const char* src) { return exactly<'+'>(src); };
    const char* kwd_minus(const char* src) { return exactly<'-'>(src); };
    const char* kwd_slash(const char* src) { return exactly<'/'>(src); };

    //####################################
    // implement some function that do exist in the standard
    // but those are locale aware which brought some trouble
    // this even seems to improve performance by quite a bit
    //####################################

    bool is_alpha(const char& chr)
    {
      return unsigned(chr - 'A') <= 'Z' - 'A' ||
             unsigned(chr - 'a') <= 'z' - 'a';
    }

    bool is_space(const char& chr)
    {
      // adapted the technique from is_alpha
      return chr == ' ' || unsigned(chr - '\t') <= '\r' - '\t';
    }

    bool is_digit(const char& chr)
    {
      // adapted the technique from is_alpha
      return unsigned(chr - '0') <= '9' - '0';
    }

    bool is_number(const char& chr)
    {
      // adapted the technique from is_alpha
      return is_digit(chr) || chr == '-' || chr == '+';
    }

    bool is_xdigit(const char& chr)
    {
      // adapted the technique from is_alpha
      return unsigned(chr - '0') <= '9' - '0' ||
             unsigned(chr - 'a') <= 'f' - 'a' ||
             unsigned(chr - 'A') <= 'F' - 'A';
    }

    bool is_punct(const char& chr)
    {
      // locale independent
      return chr == '.';
    }

    bool is_alnum(const char& chr)
    {
      return is_alpha(chr) || is_digit(chr);
    }

    // check if char is outside ascii range
    bool is_unicode(const char& chr)
    {
      // check for unicode range
      return unsigned(chr) > 127;
    }

    // check if char is outside ascii range
    // but with specific ranges (copied from Ruby Sass)
    bool is_nonascii(const char& chr)
    {
      unsigned int cmp = unsigned(chr);
      return (
        (cmp >= 128 && cmp <= 15572911) ||
        (cmp >= 15630464 && cmp <= 15712189) ||
        (cmp >= 4036001920)
      );
    }

    // check if char is within a reduced ascii range
    // valid in a uri (copied from Ruby Sass)
    bool is_uri_character(const char& chr)
    {
      unsigned int cmp = unsigned(chr);
      return (cmp > 41 && cmp < 127) ||
             cmp == ':' || cmp == '/';
    }

    // check if char is within a reduced ascii range
    // valid for escaping (copied from Ruby Sass)
    bool is_escapable_character(const char& chr)
    {
      unsigned int cmp = unsigned(chr);
      return cmp > 31 && cmp < 127;
    }

    // Match word character (look ahead)
    bool is_character(const char& chr)
    {
      // valid alpha, numeric or unicode char (plus hyphen)
      return is_alnum(chr) || is_unicode(chr) || chr == '-';
    }

    //####################################
    // BASIC CLASS MATCHERS
    //####################################

    // create matchers that advance the position
    const char* space(const char* src) { return is_space(*src) ? src + 1 : 0; }
    const char* alpha(const char* src) { return is_alpha(*src) ? src + 1 : 0; }
    const char* unicode(const char* src) { return is_unicode(*src) ? src + 1 : 0; }
    const char* nonascii(const char* src) { return is_nonascii(*src) ? src + 1 : 0; }
    const char* digit(const char* src) { return is_digit(*src) ? src + 1 : 0; }
    const char* xdigit(const char* src) { return is_xdigit(*src) ? src + 1 : 0; }
    const char* alnum(const char* src) { return is_alnum(*src) ? src + 1 : 0; }
    const char* punct(const char* src) { return is_punct(*src) ? src + 1 : 0; }
    const char* hyphen(const char* src) { return *src && *src == '-' ? src + 1 : 0; }
    const char* character(const char* src) { return is_character(*src) ? src + 1 : 0; }
    const char* uri_character(const char* src) { return is_uri_character(*src) ? src + 1 : 0; }
    const char* escapable_character(const char* src) { return is_escapable_character(*src) ? src + 1 : 0; }

    // Match multiple ctype characters.
    const char* spaces(const char* src) { return one_plus<space>(src); }
    const char* digits(const char* src) { return one_plus<digit>(src); }
    const char* hyphens(const char* src) { return one_plus<hyphen>(src); }

    // Whitespace handling.
    const char* no_spaces(const char* src) { return negate< space >(src); }
    const char* optional_spaces(const char* src) { return zero_plus< space >(src); }

    // Match any single character.
    const char* any_char(const char* src) { return *src ? src + 1 : src; }

    // Match word boundary (zero-width lookahead).
    const char* word_boundary(const char* src) { return is_character(*src) || *src == '#' ? 0 : src; }

    // Match linefeed /(?:\n|\r\n?)/
    const char* re_linebreak(const char* src)
    {
      // end of file or unix linefeed return here
      if (*src == 0 || *src == '\n') return src + 1;
      // a carriage return may optionally be followed by a linefeed
      if (*src == '\r') return *(src + 1) == '\n' ? src + 2 : src + 1;
      // no linefeed
      return 0;
    }

    // Assert string boundaries (/\Z|\z|\A/)
    // This is a zero-width positive lookahead
    const char* end_of_line(const char* src)
    {
      // end of file or unix linefeed return here
      return *src == 0 || *src == '\n' || *src == '\r' ? src : 0;
    }

    // Assert end_of_file boundary (/\z/)
    // This is a zero-width positive lookahead
    const char* end_of_file(const char* src)
    {
      // end of file or unix linefeed return here
      return *src == 0 ? src : 0;
    }

  }
}

namespace Sass {

  const std::string traces_to_string(Backtraces traces, std::string indent) {

    std::stringstream ss;
    std::string cwd(File::get_cwd());

    bool first = true;
    size_t i_beg = traces.size() - 1;
    size_t i_end = std::string::npos;
    for (size_t i = i_beg; i != i_end; i --) {

      const Backtrace& trace = traces[i];

      // make path relative to the current directory
      std::string rel_path(File::abs2rel(trace.pstate.path, cwd, cwd));

      // skip functions on error cases (unsure why ruby sass does this)
      // if (trace.caller.substr(0, 6) == ", in f") continue;

      if (first) {
        ss << indent;
        ss << "on line ";
        ss << trace.pstate.line + 1;
        ss << ":";
        ss << trace.pstate.column + 1;
        ss << " of " << rel_path;
        // ss << trace.caller;
        first = false;
      } else {
        ss << trace.caller;
        ss << std::endl;
        ss << indent;
        ss << "from line ";
        ss << trace.pstate.line + 1;
        ss << ":";
        ss << trace.pstate.column + 1;
        ss << " of " << rel_path;
      }

    }

    ss << std::endl;
    return ss.str();

  }

};
