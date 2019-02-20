#pragma once
#ifndef SASS_MAPPING_H
#define SASS_MAPPING_H

#ifndef SASS_POSITION_H
#define SASS_POSITION_H

#include <string>
#include <cstring>
// #include <iostream>

namespace Sass {


  class Offset {

    public: // c-tor
      Offset(const char chr);
      Offset(const char* string);
      Offset(const std::string& text);
      Offset(const size_t line, const size_t column);

      // return new position, incremented by the given string
      Offset add(const char* begin, const char* end);
      Offset inc(const char* begin, const char* end) const;

      // init/create instance from const char substring
      static Offset init(const char* beg, const char* end);

    public: // overload operators for position
      void operator+= (const Offset &pos);
      bool operator== (const Offset &pos) const;
      bool operator!= (const Offset &pos) const;
      Offset operator+ (const Offset &off) const;
      Offset operator- (const Offset &off) const;

    public: // overload output stream operator
      // friend std::ostream& operator<<(std::ostream& strm, const Offset& off);

    public:
      Offset off() { return *this; }

    public:
      size_t line;
      size_t column;

  };

  class Position : public Offset {

    public: // c-tor
      Position(const size_t file); // line(0), column(0)
      Position(const size_t file, const Offset& offset);
      Position(const size_t line, const size_t column); // file(-1)
      Position(const size_t file, const size_t line, const size_t column);

    public: // overload operators for position
      void operator+= (const Offset &off);
      bool operator== (const Position &pos) const;
      bool operator!= (const Position &pos) const;
      const Position operator+ (const Offset &off) const;
      const Offset operator- (const Offset &off) const;
      // return new position, incremented by the given string
      Position add(const char* begin, const char* end);
      Position inc(const char* begin, const char* end) const;

    public: // overload output stream operator
      // friend std::ostream& operator<<(std::ostream& strm, const Position& pos);

    public:
      size_t file;

  };

  // Token type for representing lexed chunks of text
  class Token {
  public:
    const char* prefix;
    const char* begin;
    const char* end;

    Token()
    : prefix(0), begin(0), end(0) { }
    Token(const char* b, const char* e)
    : prefix(b), begin(b), end(e) { }
    Token(const char* str)
    : prefix(str), begin(str), end(str + strlen(str)) { }
    Token(const char* p, const char* b, const char* e)
    : prefix(p), begin(b), end(e) { }

    size_t length()    const { return end - begin; }
    std::string ws_before() const { return std::string(prefix, begin); }
    const std::string to_string() const { return std::string(begin, end); }
    std::string time_wspace() const {
      std::string str(to_string());
      std::string whitespaces(" \t\f\v\n\r");
      return str.erase(str.find_last_not_of(whitespaces)+1);
    }

    operator bool()        { return begin && end && begin >= end; }
    operator std::string() { return to_string(); }

    bool operator==(Token t)  { return to_string() == t.to_string(); }
  };

  class ParserState : public Position {

    public: // c-tor
      ParserState(const char* path, const char* src = 0, const size_t file = std::string::npos);
      ParserState(const char* path, const char* src, const Position& position, Offset offset = Offset(0, 0));
      ParserState(const char* path, const char* src, const Token& token, const Position& position, Offset offset = Offset(0, 0));

    public: // down casts
      Offset off() { return *this; }
      Position pos() { return *this; }
      ParserState pstate() { return *this; }

    public:
      const char* path;
      const char* src;
      Offset offset;
      Token token;

  };

}

#endif

namespace Sass {

  struct Mapping {
    Position original_position;
    Position generated_position;

    Mapping(const Position& original_position, const Position& generated_position)
    : original_position(original_position), generated_position(generated_position) { }
  };

}

#endif
#ifndef SASS_C2AST_H
#define SASS_C2AST_H

#ifndef SASS_BACKTRACE_H
#define SASS_BACKTRACE_H

#include <vector>
#include <sstream>
#ifndef SASS_FILE_H
#define SASS_FILE_H


#ifndef SASS_C_CONTEXT_H
#define SASS_C_CONTEXT_H

#include <stddef.h>
#include <stdbool.h>
#ifndef SASS_BASE_H
#define SASS_BASE_H

// #define DEBUG
// #define DEBUG_SHARED_PTR

#ifdef _MSC_VER
  #pragma warning(disable : 4503)
  #ifndef _SCL_SECURE_NO_WARNINGS
    #define _SCL_SECURE_NO_WARNINGS
  #endif
  #ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS
  #endif
  #ifndef _CRT_NONSTDC_NO_DEPRECATE
    #define _CRT_NONSTDC_NO_DEPRECATE
  #endif
#endif


#ifdef __GNUC__
  #define DEPRECATED(func) func __attribute__ ((deprecated))
#elif defined(_MSC_VER)
  #define DEPRECATED(func) __declspec(deprecated) func
#else
  #pragma message("WARNING: You need to implement DEPRECATED for this compiler")
  #define DEPRECATED(func) func
#endif

#ifdef _WIN32

  /* You should define ADD_EXPORTS *only* when building the DLL. */
  #ifdef ADD_EXPORTS
    #define ADDAPI __declspec(dllexport)
    #define ADDCALL __cdecl
  #else
    #define ADDAPI
    #define ADDCALL
  #endif

#else /* _WIN32 not defined. */

  /* Define with no value on non-Windows OSes. */
  #define ADDAPI
  #define ADDCALL

#endif

/* Make sure functions are exported with C linkage under C++ compilers. */
#ifdef __cplusplus
extern "C" {
#endif


// Different render styles
enum Sass_Output_Style {
  SASS_STYLE_NESTED,
  SASS_STYLE_EXPANDED,
  SASS_STYLE_COMPACT,
  SASS_STYLE_COMPRESSED,
  // only used internaly
  SASS_STYLE_INSPECT,
  SASS_STYLE_TO_SASS
};

// to allocate buffer to be filled
ADDAPI void* ADDCALL sass_alloc_memory(size_t size);
// to allocate a buffer from existing string
ADDAPI char* ADDCALL sass_copy_c_string(const char* str);
// to free overtaken memory when done
ADDAPI void ADDCALL sass_free_memory(void* ptr);

// Some convenient string helper function
ADDAPI char* ADDCALL sass_string_quote (const char* str, const char quote_mark);
ADDAPI char* ADDCALL sass_string_unquote (const char* str);

// Implemented sass language version
// Hardcoded version 3.4 for time being
ADDAPI const char* ADDCALL libsass_version(void);

// Get compiled libsass language
ADDAPI const char* ADDCALL libsass_language_version(void);

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif
#ifndef SASS_C_VALUES_H
#define SASS_C_VALUES_H


#ifdef __cplusplus
extern "C" {
#endif


// Forward declaration
union Sass_Value;

// Type for Sass values
enum Sass_Tag {
  SASS_BOOLEAN,
  SASS_NUMBER,
  SASS_COLOR,
  SASS_STRING,
  SASS_LIST,
  SASS_MAP,
  SASS_NULL,
  SASS_ERROR,
  SASS_WARNING
};

// Tags for denoting Sass list separators
enum Sass_Separator {
  SASS_COMMA,
  SASS_SPACE,
  // only used internally to represent a hash map before evaluation
  // otherwise we would be too early to check for duplicate keys
  SASS_HASH
};

// Value Operators
enum Sass_OP {
  AND, OR,                   // logical connectives
  EQ, NEQ, GT, GTE, LT, LTE, // arithmetic relations
  ADD, SUB, MUL, DIV, MOD,   // arithmetic functions
  NUM_OPS                    // so we know how big to make the op table
};

// Creator functions for all value types
ADDAPI union Sass_Value* ADDCALL sass_make_null    (void);
ADDAPI union Sass_Value* ADDCALL sass_make_boolean (bool val);
ADDAPI union Sass_Value* ADDCALL sass_make_string  (const char* val);
ADDAPI union Sass_Value* ADDCALL sass_make_qstring (const char* val);
ADDAPI union Sass_Value* ADDCALL sass_make_number  (double val, const char* unit);
ADDAPI union Sass_Value* ADDCALL sass_make_color   (double r, double g, double b, double a);
ADDAPI union Sass_Value* ADDCALL sass_make_list    (size_t len, enum Sass_Separator sep, bool is_bracketed);
ADDAPI union Sass_Value* ADDCALL sass_make_map     (size_t len);
ADDAPI union Sass_Value* ADDCALL sass_make_error   (const char* msg);
ADDAPI union Sass_Value* ADDCALL sass_make_warning (const char* msg);

// Generic destructor function for all types
// Will release memory of all associated Sass_Values
// Means we will delete recursively for lists and maps
ADDAPI void ADDCALL sass_delete_value (union Sass_Value* val);

// Make a deep cloned copy of the given sass value
ADDAPI union Sass_Value* ADDCALL sass_clone_value (const union Sass_Value* val);

// Execute an operation for two Sass_Values and return the result as a Sass_Value too
ADDAPI union Sass_Value* ADDCALL sass_value_op (enum Sass_OP op, const union Sass_Value* a, const union Sass_Value* b);

// Stringify a Sass_Values and also return the result as a Sass_Value (of type STRING)
ADDAPI union Sass_Value* ADDCALL sass_value_stringify (const union Sass_Value* a, bool compressed, int precision);

// Return the sass tag for a generic sass value
// Check is needed before accessing specific values!
ADDAPI enum Sass_Tag ADDCALL sass_value_get_tag (const union Sass_Value* v);

// Check value to be of a specific type
// Can also be used before accessing properties!
ADDAPI bool ADDCALL sass_value_is_null (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_number (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_string (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_boolean (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_color (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_list (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_map (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_error (const union Sass_Value* v);
ADDAPI bool ADDCALL sass_value_is_warning (const union Sass_Value* v);

// Getters and setters for Sass_Number
ADDAPI double ADDCALL sass_number_get_value (const union Sass_Value* v);
ADDAPI void ADDCALL sass_number_set_value (union Sass_Value* v, double value);
ADDAPI const char* ADDCALL sass_number_get_unit (const union Sass_Value* v);
ADDAPI void ADDCALL sass_number_set_unit (union Sass_Value* v, char* unit);

// Getters and setters for Sass_String
ADDAPI const char* ADDCALL sass_string_get_value (const union Sass_Value* v);
ADDAPI void ADDCALL sass_string_set_value (union Sass_Value* v, char* value);
ADDAPI bool ADDCALL sass_string_is_quoted(const union Sass_Value* v);
ADDAPI void ADDCALL sass_string_set_quoted(union Sass_Value* v, bool quoted);

// Getters and setters for Sass_Boolean
ADDAPI bool ADDCALL sass_boolean_get_value (const union Sass_Value* v);
ADDAPI void ADDCALL sass_boolean_set_value (union Sass_Value* v, bool value);

// Getters and setters for Sass_Color
ADDAPI double ADDCALL sass_color_get_r (const union Sass_Value* v);
ADDAPI void ADDCALL sass_color_set_r (union Sass_Value* v, double r);
ADDAPI double ADDCALL sass_color_get_g (const union Sass_Value* v);
ADDAPI void ADDCALL sass_color_set_g (union Sass_Value* v, double g);
ADDAPI double ADDCALL sass_color_get_b (const union Sass_Value* v);
ADDAPI void ADDCALL sass_color_set_b (union Sass_Value* v, double b);
ADDAPI double ADDCALL sass_color_get_a (const union Sass_Value* v);
ADDAPI void ADDCALL sass_color_set_a (union Sass_Value* v, double a);

// Getter for the number of items in list
ADDAPI size_t ADDCALL sass_list_get_length (const union Sass_Value* v);
// Getters and setters for Sass_List
ADDAPI enum Sass_Separator ADDCALL sass_list_get_separator (const union Sass_Value* v);
ADDAPI void ADDCALL sass_list_set_separator (union Sass_Value* v, enum Sass_Separator value);
ADDAPI bool ADDCALL sass_list_get_is_bracketed (const union Sass_Value* v);
ADDAPI void ADDCALL sass_list_set_is_bracketed (union Sass_Value* v, bool value);
// Getters and setters for Sass_List values
ADDAPI union Sass_Value* ADDCALL sass_list_get_value (const union Sass_Value* v, size_t i);
ADDAPI void ADDCALL sass_list_set_value (union Sass_Value* v, size_t i, union Sass_Value* value);

// Getter for the number of items in map
ADDAPI size_t ADDCALL sass_map_get_length (const union Sass_Value* v);
// Getters and setters for Sass_Map keys and values
ADDAPI union Sass_Value* ADDCALL sass_map_get_key (const union Sass_Value* v, size_t i);
ADDAPI void ADDCALL sass_map_set_key (union Sass_Value* v, size_t i, union Sass_Value*);
ADDAPI union Sass_Value* ADDCALL sass_map_get_value (const union Sass_Value* v, size_t i);
ADDAPI void ADDCALL sass_map_set_value (union Sass_Value* v, size_t i, union Sass_Value*);

// Getters and setters for Sass_Error
ADDAPI char* ADDCALL sass_error_get_message (const union Sass_Value* v);
ADDAPI void ADDCALL sass_error_set_message (union Sass_Value* v, char* msg);

// Getters and setters for Sass_Warning
ADDAPI char* ADDCALL sass_warning_get_message (const union Sass_Value* v);
ADDAPI void ADDCALL sass_warning_set_message (union Sass_Value* v, char* msg);

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif
#ifndef SASS_C_FUNCTIONS_H
#define SASS_C_FUNCTIONS_H


#ifdef __cplusplus
extern "C" {
#endif


// Forward declaration
struct Sass_Env;
struct Sass_Callee;
struct Sass_Import;
struct Sass_Options;
struct Sass_Compiler;
struct Sass_Importer;
struct Sass_Function;

// Typedef helpers for callee lists
typedef struct Sass_Env (*Sass_Env_Frame);
// Typedef helpers for callee lists
typedef struct Sass_Callee (*Sass_Callee_Entry);
// Typedef helpers for import lists
typedef struct Sass_Import (*Sass_Import_Entry);
typedef struct Sass_Import* (*Sass_Import_List);
// Typedef helpers for custom importer lists
typedef struct Sass_Importer (*Sass_Importer_Entry);
typedef struct Sass_Importer* (*Sass_Importer_List);
// Typedef defining importer signature and return type
typedef Sass_Import_List (*Sass_Importer_Fn)
  (const char* url, Sass_Importer_Entry cb, struct Sass_Compiler* compiler);

// Typedef helpers for custom functions lists
typedef struct Sass_Function (*Sass_Function_Entry);
typedef struct Sass_Function* (*Sass_Function_List);
// Typedef defining function signature and return type
typedef union Sass_Value* (*Sass_Function_Fn)
  (const union Sass_Value*, Sass_Function_Entry cb, struct Sass_Compiler* compiler);

// Type of function calls
enum Sass_Callee_Type {
  SASS_CALLEE_MIXIN,
  SASS_CALLEE_FUNCTION,
  SASS_CALLEE_C_FUNCTION,
};

// Creator for sass custom importer return argument list
ADDAPI Sass_Importer_List ADDCALL sass_make_importer_list (size_t length);
ADDAPI Sass_Importer_Entry ADDCALL sass_importer_get_list_entry (Sass_Importer_List list, size_t idx);
ADDAPI void ADDCALL sass_importer_set_list_entry (Sass_Importer_List list, size_t idx, Sass_Importer_Entry entry);
ADDAPI void ADDCALL sass_delete_importer_list (Sass_Importer_List list);


// Creators for custom importer callback (with some additional pointer)
// The pointer is mostly used to store the callback into the actual binding
ADDAPI Sass_Importer_Entry ADDCALL sass_make_importer (Sass_Importer_Fn importer, double priority, void* cookie);

// Getters for import function descriptors
ADDAPI Sass_Importer_Fn ADDCALL sass_importer_get_function (Sass_Importer_Entry cb);
ADDAPI double ADDCALL sass_importer_get_priority (Sass_Importer_Entry cb);
ADDAPI void* ADDCALL sass_importer_get_cookie (Sass_Importer_Entry cb);

// Deallocator for associated memory
ADDAPI void ADDCALL sass_delete_importer (Sass_Importer_Entry cb);

// Creator for sass custom importer return argument list
ADDAPI Sass_Import_List ADDCALL sass_make_import_list (size_t length);
// Creator for a single import entry returned by the custom importer inside the list
ADDAPI Sass_Import_Entry ADDCALL sass_make_import_entry (const char* path, char* source, char* srcmap);
ADDAPI Sass_Import_Entry ADDCALL sass_make_import (const char* imp_path, const char* abs_base, char* source, char* srcmap);
// set error message to abort import and to print out a message (path from existing object is used in output)
ADDAPI Sass_Import_Entry ADDCALL sass_import_set_error(Sass_Import_Entry import, const char* message, size_t line, size_t col);

// Setters to insert an entry into the import list (you may also use [] access directly)
// Since we are dealing with pointers they should have a guaranteed and fixed size
ADDAPI void ADDCALL sass_import_set_list_entry (Sass_Import_List list, size_t idx, Sass_Import_Entry entry);
ADDAPI Sass_Import_Entry ADDCALL sass_import_get_list_entry (Sass_Import_List list, size_t idx);

// Getters for callee entry
ADDAPI const char* ADDCALL sass_callee_get_name (Sass_Callee_Entry);
ADDAPI const char* ADDCALL sass_callee_get_path (Sass_Callee_Entry);
ADDAPI size_t ADDCALL sass_callee_get_line (Sass_Callee_Entry);
ADDAPI size_t ADDCALL sass_callee_get_column (Sass_Callee_Entry);
ADDAPI enum Sass_Callee_Type ADDCALL sass_callee_get_type (Sass_Callee_Entry);
ADDAPI Sass_Env_Frame ADDCALL sass_callee_get_env (Sass_Callee_Entry);

// Getters and Setters for environments (lexical, local and global)
ADDAPI union Sass_Value* ADDCALL sass_env_get_lexical (Sass_Env_Frame, const char*);
ADDAPI void ADDCALL sass_env_set_lexical (Sass_Env_Frame, const char*, union Sass_Value*);
ADDAPI union Sass_Value* ADDCALL sass_env_get_local (Sass_Env_Frame, const char*);
ADDAPI void ADDCALL sass_env_set_local (Sass_Env_Frame, const char*, union Sass_Value*);
ADDAPI union Sass_Value* ADDCALL sass_env_get_global (Sass_Env_Frame, const char*);
ADDAPI void ADDCALL sass_env_set_global (Sass_Env_Frame, const char*, union Sass_Value*);

// Getters for import entry
ADDAPI const char* ADDCALL sass_import_get_imp_path (Sass_Import_Entry);
ADDAPI const char* ADDCALL sass_import_get_abs_path (Sass_Import_Entry);
ADDAPI const char* ADDCALL sass_import_get_source (Sass_Import_Entry);
ADDAPI const char* ADDCALL sass_import_get_srcmap (Sass_Import_Entry);
// Explicit functions to take ownership of these items
// The property on our struct will be reset to NULL
ADDAPI char* ADDCALL sass_import_take_source (Sass_Import_Entry);
ADDAPI char* ADDCALL sass_import_take_srcmap (Sass_Import_Entry);
// Getters from import error entry
ADDAPI size_t ADDCALL sass_import_get_error_line (Sass_Import_Entry);
ADDAPI size_t ADDCALL sass_import_get_error_column (Sass_Import_Entry);
ADDAPI const char* ADDCALL sass_import_get_error_message (Sass_Import_Entry);

// Deallocator for associated memory (incl. entries)
ADDAPI void ADDCALL sass_delete_import_list (Sass_Import_List);
// Just in case we have some stray import structs
ADDAPI void ADDCALL sass_delete_import (Sass_Import_Entry);



// Creators for sass function list and function descriptors
ADDAPI Sass_Function_List ADDCALL sass_make_function_list (size_t length);
ADDAPI Sass_Function_Entry ADDCALL sass_make_function (const char* signature, Sass_Function_Fn cb, void* cookie);
ADDAPI void ADDCALL sass_delete_function (Sass_Function_Entry entry);
ADDAPI void ADDCALL sass_delete_function_list (Sass_Function_List list);

// Setters and getters for callbacks on function lists
ADDAPI Sass_Function_Entry ADDCALL sass_function_get_list_entry(Sass_Function_List list, size_t pos);
ADDAPI void ADDCALL sass_function_set_list_entry(Sass_Function_List list, size_t pos, Sass_Function_Entry cb);

// Getters for custom function descriptors
ADDAPI const char* ADDCALL sass_function_get_signature (Sass_Function_Entry cb);
ADDAPI Sass_Function_Fn ADDCALL sass_function_get_function (Sass_Function_Entry cb);
ADDAPI void* ADDCALL sass_function_get_cookie (Sass_Function_Entry cb);


#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif

#ifdef __cplusplus
extern "C" {
#endif


// Forward declaration
struct Sass_Compiler;

// Forward declaration
struct Sass_Options; // base struct
struct Sass_Context; // : Sass_Options
struct Sass_File_Context; // : Sass_Context
struct Sass_Data_Context; // : Sass_Context

// Compiler states
enum Sass_Compiler_State {
  SASS_COMPILER_CREATED,
  SASS_COMPILER_PARSED,
  SASS_COMPILER_EXECUTED
};

// Create and initialize an option struct
ADDAPI struct Sass_Options* ADDCALL sass_make_options (void);
// Create and initialize a specific context
ADDAPI struct Sass_File_Context* ADDCALL sass_make_file_context (const char* input_path);
ADDAPI struct Sass_Data_Context* ADDCALL sass_make_data_context (char* source_string);

// Call the compilation step for the specific context
ADDAPI int ADDCALL sass_compile_file_context (struct Sass_File_Context* ctx);
ADDAPI int ADDCALL sass_compile_data_context (struct Sass_Data_Context* ctx);

// Create a sass compiler instance for more control
ADDAPI struct Sass_Compiler* ADDCALL sass_make_file_compiler (struct Sass_File_Context* file_ctx);
ADDAPI struct Sass_Compiler* ADDCALL sass_make_data_compiler (struct Sass_Data_Context* data_ctx);

// Execute the different compilation steps individually
// Useful if you only want to query the included files
ADDAPI int ADDCALL sass_compiler_parse(struct Sass_Compiler* compiler);
ADDAPI int ADDCALL sass_compiler_execute(struct Sass_Compiler* compiler);

// Release all memory allocated with the compiler
// This does _not_ include any contexts or options
ADDAPI void ADDCALL sass_delete_compiler(struct Sass_Compiler* compiler);
ADDAPI void ADDCALL sass_delete_options(struct Sass_Options* options);

// Release all memory allocated and also ourself
ADDAPI void ADDCALL sass_delete_file_context (struct Sass_File_Context* ctx);
ADDAPI void ADDCALL sass_delete_data_context (struct Sass_Data_Context* ctx);

// Getters for context from specific implementation
ADDAPI struct Sass_Context* ADDCALL sass_file_context_get_context (struct Sass_File_Context* file_ctx);
ADDAPI struct Sass_Context* ADDCALL sass_data_context_get_context (struct Sass_Data_Context* data_ctx);

// Getters for Context_Options from Sass_Context
ADDAPI struct Sass_Options* ADDCALL sass_context_get_options (struct Sass_Context* ctx);
ADDAPI struct Sass_Options* ADDCALL sass_file_context_get_options (struct Sass_File_Context* file_ctx);
ADDAPI struct Sass_Options* ADDCALL sass_data_context_get_options (struct Sass_Data_Context* data_ctx);
ADDAPI void ADDCALL sass_file_context_set_options (struct Sass_File_Context* file_ctx, struct Sass_Options* opt);
ADDAPI void ADDCALL sass_data_context_set_options (struct Sass_Data_Context* data_ctx, struct Sass_Options* opt);


// Getters for Context_Option values
ADDAPI int ADDCALL sass_option_get_precision (struct Sass_Options* options);
ADDAPI enum Sass_Output_Style ADDCALL sass_option_get_output_style (struct Sass_Options* options);
ADDAPI bool ADDCALL sass_option_get_source_comments (struct Sass_Options* options);
ADDAPI bool ADDCALL sass_option_get_source_map_embed (struct Sass_Options* options);
ADDAPI bool ADDCALL sass_option_get_source_map_contents (struct Sass_Options* options);
ADDAPI bool ADDCALL sass_option_get_source_map_file_urls (struct Sass_Options* options);
ADDAPI bool ADDCALL sass_option_get_omit_source_map_url (struct Sass_Options* options);
ADDAPI bool ADDCALL sass_option_get_is_indented_syntax_src (struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_indent (struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_linefeed (struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_input_path (struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_output_path (struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_source_map_file (struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_source_map_root (struct Sass_Options* options);
ADDAPI Sass_Importer_List ADDCALL sass_option_get_c_headers (struct Sass_Options* options);
ADDAPI Sass_Importer_List ADDCALL sass_option_get_c_importers (struct Sass_Options* options);
ADDAPI Sass_Function_List ADDCALL sass_option_get_c_functions (struct Sass_Options* options);

// Setters for Context_Option values
ADDAPI void ADDCALL sass_option_set_precision (struct Sass_Options* options, int precision);
ADDAPI void ADDCALL sass_option_set_output_style (struct Sass_Options* options, enum Sass_Output_Style output_style);
ADDAPI void ADDCALL sass_option_set_source_comments (struct Sass_Options* options, bool source_comments);
ADDAPI void ADDCALL sass_option_set_source_map_embed (struct Sass_Options* options, bool source_map_embed);
ADDAPI void ADDCALL sass_option_set_source_map_contents (struct Sass_Options* options, bool source_map_contents);
ADDAPI void ADDCALL sass_option_set_source_map_file_urls (struct Sass_Options* options, bool source_map_file_urls);
ADDAPI void ADDCALL sass_option_set_omit_source_map_url (struct Sass_Options* options, bool omit_source_map_url);
ADDAPI void ADDCALL sass_option_set_is_indented_syntax_src (struct Sass_Options* options, bool is_indented_syntax_src);
ADDAPI void ADDCALL sass_option_set_indent (struct Sass_Options* options, const char* indent);
ADDAPI void ADDCALL sass_option_set_linefeed (struct Sass_Options* options, const char* linefeed);
ADDAPI void ADDCALL sass_option_set_input_path (struct Sass_Options* options, const char* input_path);
ADDAPI void ADDCALL sass_option_set_output_path (struct Sass_Options* options, const char* output_path);
ADDAPI void ADDCALL sass_option_set_plugin_path (struct Sass_Options* options, const char* plugin_path);
ADDAPI void ADDCALL sass_option_set_include_path (struct Sass_Options* options, const char* include_path);
ADDAPI void ADDCALL sass_option_set_source_map_file (struct Sass_Options* options, const char* source_map_file);
ADDAPI void ADDCALL sass_option_set_source_map_root (struct Sass_Options* options, const char* source_map_root);
ADDAPI void ADDCALL sass_option_set_c_headers (struct Sass_Options* options, Sass_Importer_List c_headers);
ADDAPI void ADDCALL sass_option_set_c_importers (struct Sass_Options* options, Sass_Importer_List c_importers);
ADDAPI void ADDCALL sass_option_set_c_functions (struct Sass_Options* options, Sass_Function_List c_functions);


// Getters for Sass_Context values
ADDAPI const char* ADDCALL sass_context_get_output_string (struct Sass_Context* ctx);
ADDAPI int ADDCALL sass_context_get_error_status (struct Sass_Context* ctx);
ADDAPI const char* ADDCALL sass_context_get_error_json (struct Sass_Context* ctx);
ADDAPI const char* ADDCALL sass_context_get_error_text (struct Sass_Context* ctx);
ADDAPI const char* ADDCALL sass_context_get_error_message (struct Sass_Context* ctx);
ADDAPI const char* ADDCALL sass_context_get_error_file (struct Sass_Context* ctx);
ADDAPI const char* ADDCALL sass_context_get_error_src (struct Sass_Context* ctx);
ADDAPI size_t ADDCALL sass_context_get_error_line (struct Sass_Context* ctx);
ADDAPI size_t ADDCALL sass_context_get_error_column (struct Sass_Context* ctx);
ADDAPI const char* ADDCALL sass_context_get_source_map_string (struct Sass_Context* ctx);
ADDAPI char** ADDCALL sass_context_get_included_files (struct Sass_Context* ctx);

// Getters for options include path array
ADDAPI size_t ADDCALL sass_option_get_include_path_size(struct Sass_Options* options);
ADDAPI const char* ADDCALL sass_option_get_include_path(struct Sass_Options* options, size_t i);

// Calculate the size of the stored null terminated array
ADDAPI size_t ADDCALL sass_context_get_included_files_size (struct Sass_Context* ctx);

// Take ownership of memory (value on context is set to 0)
ADDAPI char* ADDCALL sass_context_take_error_json (struct Sass_Context* ctx);
ADDAPI char* ADDCALL sass_context_take_error_text (struct Sass_Context* ctx);
ADDAPI char* ADDCALL sass_context_take_error_message (struct Sass_Context* ctx);
ADDAPI char* ADDCALL sass_context_take_error_file (struct Sass_Context* ctx);
ADDAPI char* ADDCALL sass_context_take_output_string (struct Sass_Context* ctx);
ADDAPI char* ADDCALL sass_context_take_source_map_string (struct Sass_Context* ctx);
ADDAPI char** ADDCALL sass_context_take_included_files (struct Sass_Context* ctx);

// Getters for Sass_Compiler options
ADDAPI enum Sass_Compiler_State ADDCALL sass_compiler_get_state(struct Sass_Compiler* compiler);
ADDAPI struct Sass_Context* ADDCALL sass_compiler_get_context(struct Sass_Compiler* compiler);
ADDAPI struct Sass_Options* ADDCALL sass_compiler_get_options(struct Sass_Compiler* compiler);
ADDAPI size_t ADDCALL sass_compiler_get_import_stack_size(struct Sass_Compiler* compiler);
ADDAPI Sass_Import_Entry ADDCALL sass_compiler_get_last_import(struct Sass_Compiler* compiler);
ADDAPI Sass_Import_Entry ADDCALL sass_compiler_get_import_entry(struct Sass_Compiler* compiler, size_t idx);
ADDAPI size_t ADDCALL sass_compiler_get_callee_stack_size(struct Sass_Compiler* compiler);
ADDAPI Sass_Callee_Entry ADDCALL sass_compiler_get_last_callee(struct Sass_Compiler* compiler);
ADDAPI Sass_Callee_Entry ADDCALL sass_compiler_get_callee_entry(struct Sass_Compiler* compiler, size_t idx);

// Push function for paths (no manipulation support for now)
ADDAPI void ADDCALL sass_option_push_plugin_path (struct Sass_Options* options, const char* path);
ADDAPI void ADDCALL sass_option_push_include_path (struct Sass_Options* options, const char* path);

// Resolve a file via the given include paths in the sass option struct
// find_file looks for the exact file name while find_include does a regular sass include
ADDAPI char* ADDCALL sass_find_file (const char* path, struct Sass_Options* opt);
ADDAPI char* ADDCALL sass_find_include (const char* path, struct Sass_Options* opt);

// Resolve a file relative to last import or include paths in the sass option struct
// find_file looks for the exact file name while find_include does a regular sass include
ADDAPI char* ADDCALL sass_compiler_find_file (const char* path, struct Sass_Compiler* compiler);
ADDAPI char* ADDCALL sass_compiler_find_include (const char* path, struct Sass_Compiler* compiler);

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif
#ifndef SASS_AST_FWD_DECL_H
#define SASS_AST_FWD_DECL_H

#include <map>
#include <set>
#include <deque>
#include <typeinfo>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#ifndef SASS_MEMORY_SHARED_PTR_H
#define SASS_MEMORY_SHARED_PTR_H



namespace Sass {

  class SharedPtr;

  ///////////////////////////////////////////////////////////////////////////////
  // Use macros for the allocation task, since overloading operator `new`
  // has been proven to be flaky under certain compilers (see comment below).
  ///////////////////////////////////////////////////////////////////////////////

  #ifdef DEBUG_SHARED_PTR

    #define SASS_MEMORY_NEW(Class, ...) \
      ((Class*)(new Class(__VA_ARGS__))->trace(__FILE__, __LINE__)) \

    #define SASS_MEMORY_COPY(obj) \
      ((obj)->copy(__FILE__, __LINE__)) \

    #define SASS_MEMORY_CLONE(obj) \
      ((obj)->clone(__FILE__, __LINE__)) \

  #else

    #define SASS_MEMORY_NEW(Class, ...) \
      new Class(__VA_ARGS__) \

    #define SASS_MEMORY_COPY(obj) \
      ((obj)->copy()) \

    #define SASS_MEMORY_CLONE(obj) \
      ((obj)->clone()) \

  #endif

  class SharedObj {
   public:
    SharedObj() : refcount(0), detached(false) {
      #ifdef DEBUG_SHARED_PTR
      if (taint) all.push_back(this);
      #endif
    }
    virtual ~SharedObj() {
      #ifdef DEBUG_SHARED_PTR
      all.clear();
      #endif
    }

    #ifdef DEBUG_SHARED_PTR
    static void dumpMemLeaks();
    SharedObj* trace(std::string file, size_t line) {
      this->file = file;
      this->line = line;
      return this;
    }
    std::string getDbgFile() { return file; }
    size_t getDbgLine() { return line; }
    void setDbg(bool dbg) { this->dbg = dbg; }
    size_t getRefCount() const { return refcount; }
    #endif

    static void setTaint(bool val) { taint = val; }

    virtual const std::string to_string() const = 0;
   protected:
    friend class SharedPtr;
    friend class Memory_Manager;
    size_t refcount;
    bool detached;
    static bool taint;
    #ifdef DEBUG_SHARED_PTR
    std::string file;
    size_t line;
    bool dbg = false;
    static std::vector<SharedObj*> all;
    #endif
  };

  class SharedPtr {
   public:
    SharedPtr() : node(nullptr) {}
    SharedPtr(SharedObj* ptr) : node(ptr) {
      incRefCount();
    }
    SharedPtr(const SharedPtr& obj) : SharedPtr(obj.node) {}
    ~SharedPtr() {
      decRefCount();
    }

    SharedPtr& operator=(SharedObj* other_node) {
      if (node != other_node) {
        decRefCount();
        node = other_node;
        incRefCount();
      } else if (node != nullptr) {
        node->detached = false;
      }
      return *this;
    }

    SharedPtr& operator=(const SharedPtr& obj) {
      return *this = obj.node;
    }

    // Prevents all SharedPtrs from freeing this node until it is assigned to another SharedPtr.
    SharedObj* detach() {
      if (node != nullptr) node->detached = true;
      return node;
    }

    SharedObj* obj() const { return node; }
    SharedObj* operator->() const { return node; }
    bool isNull() const { return node == nullptr; }
    operator bool() const { return node != nullptr; }

   protected:
    SharedObj* node;
    void decRefCount() {
      if (node == nullptr) return;
      --node->refcount;
      #ifdef DEBUG_SHARED_PTR
      if (node->dbg) std::cerr << "- " << node << " X " << node->refcount << " (" << this << ") " << "\n";
      #endif
      if (node->refcount == 0 && !node->detached) {
        #ifdef DEBUG_SHARED_PTR
        if (node->dbg) std::cerr << "DELETE NODE " << node << "\n";
        #endif
        delete node;
      }
    }
    void incRefCount() {
      if (node == nullptr) return;
      node->detached = false;
      ++node->refcount;
      #ifdef DEBUG_SHARED_PTR
      if (node->dbg) std::cerr << "+ " << node << " X " << node->refcount << " (" << this << ") " << "\n";
      #endif
    }
  };

  template <class T>
  class SharedImpl : private SharedPtr {
   public:
    SharedImpl() : SharedPtr(nullptr) {}

    template <class U>
    SharedImpl(U* node) :
      SharedPtr(static_cast<T*>(node)) {}

    template <class U>
    SharedImpl(const SharedImpl<U>& impl) :
      SharedImpl(impl.ptr()) {}

    template <class U>
    SharedImpl<T>& operator=(U *rhs) {
      return static_cast<SharedImpl<T>&>(
        SharedPtr::operator=(static_cast<T*>(rhs)));
    }

    template <class U>
    SharedImpl<T>& operator=(const SharedImpl<U>& rhs) {
      return static_cast<SharedImpl<T>&>(
        SharedPtr::operator=(static_cast<const SharedImpl<T>&>(rhs)));
    }

    operator const std::string() const {
      if (node) return node->to_string();
      return "[nullptr]";
    }

    using SharedPtr::isNull;
    using SharedPtr::operator bool;
    operator T*() const { return static_cast<T*>(this->obj()); }
    operator T&() const { return *static_cast<T*>(this->obj()); }
    T& operator* () const { return *static_cast<T*>(this->obj()); };
    T* operator-> () const { return static_cast<T*>(this->obj()); };
    T* ptr () const { return static_cast<T*>(this->obj()); };
    T* detach() { return static_cast<T*>(SharedPtr::detach()); }
  };

}

#endif

/////////////////////////////////////////////
// Forward declarations for the AST visitors.
/////////////////////////////////////////////
namespace Sass {

  class AST_Node;
  typedef AST_Node* AST_Node_Ptr;
  typedef AST_Node const* AST_Node_Ptr_Const;

  class Has_Block;
  typedef Has_Block* Has_Block_Ptr;
  typedef Has_Block const* Has_Block_Ptr_Const;

  class Simple_Selector;
  typedef Simple_Selector* Simple_Selector_Ptr;
  typedef Simple_Selector const* Simple_Selector_Ptr_Const;

  class Parent_Reference;
  typedef Parent_Reference* Parent_Reference_Ptr;
  typedef Parent_Reference const* Parent_Reference_Ptr_Const;

  class PreValue;
  typedef PreValue* PreValue_Ptr;
  typedef PreValue const* PreValue_Ptr_Const;
  class Block;
  typedef Block* Block_Ptr;
  typedef Block const* Block_Ptr_Const;
  class Expression;
  typedef Expression* Expression_Ptr;
  typedef Expression const* Expression_Ptr_Const;
  class Statement;
  typedef Statement* Statement_Ptr;
  typedef Statement const* Statement_Ptr_Const;
  class Value;
  typedef Value* Value_Ptr;
  typedef Value const* Value_Ptr_Const;
  class Declaration;
  typedef Declaration* Declaration_Ptr;
  typedef Declaration const* Declaration_Ptr_Const;
  class Ruleset;
  typedef Ruleset* Ruleset_Ptr;
  typedef Ruleset const* Ruleset_Ptr_Const;
  class Bubble;
  typedef Bubble* Bubble_Ptr;
  typedef Bubble const* Bubble_Ptr_Const;
  class Trace;
  typedef Trace* Trace_Ptr;
  typedef Trace const* Trace_Ptr_Const;

  class Media_Block;
  typedef Media_Block* Media_Block_Ptr;
  typedef Media_Block const* Media_Block_Ptr_Const;
  class Supports_Block;
  typedef Supports_Block* Supports_Block_Ptr;
  typedef Supports_Block const* Supports_Block_Ptr_Const;
  class Directive;
  typedef Directive* Directive_Ptr;
  typedef Directive const* Directive_Ptr_Const;


  class Keyframe_Rule;
  typedef Keyframe_Rule* Keyframe_Rule_Ptr;
  typedef Keyframe_Rule const* Keyframe_Rule_Ptr_Const;
  class At_Root_Block;
  typedef At_Root_Block* At_Root_Block_Ptr;
  typedef At_Root_Block const* At_Root_Block_Ptr_Const;
  class Assignment;
  typedef Assignment* Assignment_Ptr;
  typedef Assignment const* Assignment_Ptr_Const;

  class Import;
  typedef Import* Import_Ptr;
  typedef Import const* Import_Ptr_Const;
  class Import_Stub;
  typedef Import_Stub* Import_Stub_Ptr;
  typedef Import_Stub const* Import_Stub_Ptr_Const;
  class Warning;
  typedef Warning* Warning_Ptr;
  typedef Warning const* Warning_Ptr_Const;

  class Error;
  typedef Error* Error_Ptr;
  typedef Error const* Error_Ptr_Const;
  class Debug;
  typedef Debug* Debug_Ptr;
  typedef Debug const* Debug_Ptr_Const;
  class Comment;
  typedef Comment* Comment_Ptr;
  typedef Comment const* Comment_Ptr_Const;

  class If;
  typedef If* If_Ptr;
  typedef If const* If_Ptr_Const;
  class For;
  typedef For* For_Ptr;
  typedef For const* For_Ptr_Const;
  class Each;
  typedef Each* Each_Ptr;
  typedef Each const* Each_Ptr_Const;
  class While;
  typedef While* While_Ptr;
  typedef While const* While_Ptr_Const;
  class Return;
  typedef Return* Return_Ptr;
  typedef Return const* Return_Ptr_Const;
  class Content;
  typedef Content* Content_Ptr;
  typedef Content const* Content_Ptr_Const;
  class Extension;
  typedef Extension* Extension_Ptr;
  typedef Extension const* Extension_Ptr_Const;
  class Definition;
  typedef Definition* Definition_Ptr;
  typedef Definition const* Definition_Ptr_Const;

  class List;
  typedef List* List_Ptr;
  typedef List const* List_Ptr_Const;
  class Map;
  typedef Map* Map_Ptr;
  typedef Map const* Map_Ptr_Const;
  class Function;
  typedef Function* Function_Ptr;
  typedef Function const* Function_Ptr_Const;

  class Mixin_Call;
  typedef Mixin_Call* Mixin_Call_Ptr;
  typedef Mixin_Call const* Mixin_Call_Ptr_Const;
  class Binary_Expression;
  typedef Binary_Expression* Binary_Expression_Ptr;
  typedef Binary_Expression const* Binary_Expression_Ptr_Const;
  class Unary_Expression;
  typedef Unary_Expression* Unary_Expression_Ptr;
  typedef Unary_Expression const* Unary_Expression_Ptr_Const;
  class Function_Call;
  typedef Function_Call* Function_Call_Ptr;
  typedef Function_Call const* Function_Call_Ptr_Const;
  class Custom_Warning;
  typedef Custom_Warning* Custom_Warning_Ptr;
  typedef Custom_Warning const* Custom_Warning_Ptr_Const;
  class Custom_Error;
  typedef Custom_Error* Custom_Error_Ptr;
  typedef Custom_Error const* Custom_Error_Ptr_Const;

  class Variable;
  typedef Variable* Variable_Ptr;
  typedef Variable const* Variable_Ptr_Const;
  class Number;
  typedef Number* Number_Ptr;
  typedef Number const* Number_Ptr_Const;
  class Color;
  typedef Color* Color_Ptr;
  typedef Color const* Color_Ptr_Const;
  class Color_RGBA;
  typedef Color_RGBA* Color_RGBA_Ptr;
  typedef Color_RGBA const* Color_RGBA_Ptr_Const;
  class Color_HSLA;
  typedef Color_HSLA* Color_HSLA_Ptr;
  typedef Color_HSLA const* Color_HSLA_Ptr_Const;
  class Boolean;
  typedef Boolean* Boolean_Ptr;
  typedef Boolean const* Boolean_Ptr_Const;
  class String;
  typedef String* String_Ptr;
  typedef String const* String_Ptr_Const;

  class String_Schema;
  typedef String_Schema* String_Schema_Ptr;
  typedef String_Schema const* String_Schema_Ptr_Const;
  class String_Constant;
  typedef String_Constant* String_Constant_Ptr;
  typedef String_Constant const* String_Constant_Ptr_Const;
  class String_Quoted;
  typedef String_Quoted* String_Quoted_Ptr;
  typedef String_Quoted const* String_Quoted_Ptr_Const;

  class Media_Query;
  typedef Media_Query* Media_Query_Ptr;
  typedef Media_Query const* Media_Query_Ptr_Const;
  class Media_Query_Expression;
  typedef Media_Query_Expression* Media_Query_Expression_Ptr;
  typedef Media_Query_Expression const* Media_Query_Expression_Ptr_Const;
  class Supports_Condition;
  typedef Supports_Condition* Supports_Condition_Ptr;
  typedef Supports_Condition const* Supports_Condition_Ptr_Const;
  class Supports_Operator;
  typedef Supports_Operator* Supports_Operator_Ptr;
  typedef Supports_Operator const* Supports_Operator_Ptr_Const;
  class Supports_Negation;
  typedef Supports_Negation* Supports_Negation_Ptr;
  typedef Supports_Negation const* Supports_Negation_Ptr_Const;
  class Supports_Declaration;
  typedef Supports_Declaration* Supports_Declaration_Ptr;
  typedef Supports_Declaration const* Supports_Declaration_Ptr_Const;
  class Supports_Interpolation;
  typedef Supports_Interpolation* Supports_Interpolation_Ptr;
  typedef Supports_Interpolation const* Supports_Interpolation_Ptr_Const;


  class Null;
  typedef Null* Null_Ptr;
  typedef Null const* Null_Ptr_Const;

  class At_Root_Query;
  typedef At_Root_Query* At_Root_Query_Ptr;
  typedef At_Root_Query const* At_Root_Query_Ptr_Const;
  class Parent_Selector;
  typedef Parent_Selector* Parent_Selector_Ptr;
  typedef Parent_Selector const* Parent_Selector_Ptr_Const;
  class Parameter;
  typedef Parameter* Parameter_Ptr;
  typedef Parameter const* Parameter_Ptr_Const;
  class Parameters;
  typedef Parameters* Parameters_Ptr;
  typedef Parameters const* Parameters_Ptr_Const;
  class Argument;
  typedef Argument* Argument_Ptr;
  typedef Argument const* Argument_Ptr_Const;
  class Arguments;
  typedef Arguments* Arguments_Ptr;
  typedef Arguments const* Arguments_Ptr_Const;
  class Selector;
  typedef Selector* Selector_Ptr;
  typedef Selector const* Selector_Ptr_Const;


  class Selector_Schema;
  typedef Selector_Schema* Selector_Schema_Ptr;
  typedef Selector_Schema const* Selector_Schema_Ptr_Const;
  class Placeholder_Selector;
  typedef Placeholder_Selector* Placeholder_Selector_Ptr;
  typedef Placeholder_Selector const* Placeholder_Selector_Ptr_Const;
  class Type_Selector;
  typedef Type_Selector* Type_Selector_Ptr;
  typedef Type_Selector const* Type_Selector_Ptr_Const;
  class Class_Selector;
  typedef Class_Selector* Class_Selector_Ptr;
  typedef Class_Selector const* Class_Selector_Ptr_Const;
  class Id_Selector;
  typedef Id_Selector* Id_Selector_Ptr;
  typedef Id_Selector const* Id_Selector_Ptr_Const;
  class Attribute_Selector;
  typedef Attribute_Selector* Attribute_Selector_Ptr;
  typedef Attribute_Selector const* Attribute_Selector_Ptr_Const;

  class Pseudo_Selector;
  typedef Pseudo_Selector* Pseudo_Selector_Ptr;
  typedef Pseudo_Selector const * Pseudo_Selector_Ptr_Const;
  class Wrapped_Selector;
  typedef Wrapped_Selector* Wrapped_Selector_Ptr;
  typedef Wrapped_Selector const * Wrapped_Selector_Ptr_Const;
  class Compound_Selector;
  typedef Compound_Selector* Compound_Selector_Ptr;
  typedef Compound_Selector const * Compound_Selector_Ptr_Const;
  class Complex_Selector;
  typedef Complex_Selector* Complex_Selector_Ptr;
  typedef Complex_Selector const * Complex_Selector_Ptr_Const;
  class Selector_List;
  typedef Selector_List* Selector_List_Ptr;
  typedef Selector_List const * Selector_List_Ptr_Const;


  // common classes
  class Context;
  class Expand;
  class Eval;

  // declare classes that are instances of memory nodes
  // #define IMPL_MEM_OBJ(type) using type##_Obj = SharedImpl<type>
  #define IMPL_MEM_OBJ(type) typedef SharedImpl<type> type##_Obj

  IMPL_MEM_OBJ(AST_Node);
  IMPL_MEM_OBJ(Statement);
  IMPL_MEM_OBJ(Block);
  IMPL_MEM_OBJ(Ruleset);
  IMPL_MEM_OBJ(Bubble);
  IMPL_MEM_OBJ(Trace);
  IMPL_MEM_OBJ(Media_Block);
  IMPL_MEM_OBJ(Supports_Block);
  IMPL_MEM_OBJ(Directive);
  IMPL_MEM_OBJ(Keyframe_Rule);
  IMPL_MEM_OBJ(At_Root_Block);
  IMPL_MEM_OBJ(Declaration);
  IMPL_MEM_OBJ(Assignment);
  IMPL_MEM_OBJ(Import);
  IMPL_MEM_OBJ(Import_Stub);
  IMPL_MEM_OBJ(Warning);
  IMPL_MEM_OBJ(Error);
  IMPL_MEM_OBJ(Debug);
  IMPL_MEM_OBJ(Comment);
  IMPL_MEM_OBJ(PreValue);
  IMPL_MEM_OBJ(Has_Block);
  IMPL_MEM_OBJ(If);
  IMPL_MEM_OBJ(For);
  IMPL_MEM_OBJ(Each);
  IMPL_MEM_OBJ(While);
  IMPL_MEM_OBJ(Return);
  IMPL_MEM_OBJ(Content);
  IMPL_MEM_OBJ(Extension);
  IMPL_MEM_OBJ(Definition);
  IMPL_MEM_OBJ(Mixin_Call);
  IMPL_MEM_OBJ(Value);
  IMPL_MEM_OBJ(Expression);
  IMPL_MEM_OBJ(List);
  IMPL_MEM_OBJ(Map);
  IMPL_MEM_OBJ(Function);
  IMPL_MEM_OBJ(Binary_Expression);
  IMPL_MEM_OBJ(Unary_Expression);
  IMPL_MEM_OBJ(Function_Call);
  IMPL_MEM_OBJ(Custom_Warning);
  IMPL_MEM_OBJ(Custom_Error);
  IMPL_MEM_OBJ(Variable);
  IMPL_MEM_OBJ(Number);
  IMPL_MEM_OBJ(Color);
  IMPL_MEM_OBJ(Color_RGBA);
  IMPL_MEM_OBJ(Color_HSLA);
  IMPL_MEM_OBJ(Boolean);
  IMPL_MEM_OBJ(String_Schema);
  IMPL_MEM_OBJ(String);
  IMPL_MEM_OBJ(String_Constant);
  IMPL_MEM_OBJ(String_Quoted);
  IMPL_MEM_OBJ(Media_Query);
  IMPL_MEM_OBJ(Media_Query_Expression);
  IMPL_MEM_OBJ(Supports_Condition);
  IMPL_MEM_OBJ(Supports_Operator);
  IMPL_MEM_OBJ(Supports_Negation);
  IMPL_MEM_OBJ(Supports_Declaration);
  IMPL_MEM_OBJ(Supports_Interpolation);
  IMPL_MEM_OBJ(At_Root_Query);
  IMPL_MEM_OBJ(Null);
  IMPL_MEM_OBJ(Parent_Selector);
  IMPL_MEM_OBJ(Parent_Reference);
  IMPL_MEM_OBJ(Parameter);
  IMPL_MEM_OBJ(Parameters);
  IMPL_MEM_OBJ(Argument);
  IMPL_MEM_OBJ(Arguments);
  IMPL_MEM_OBJ(Selector);
  IMPL_MEM_OBJ(Selector_Schema);
  IMPL_MEM_OBJ(Simple_Selector);
  IMPL_MEM_OBJ(Placeholder_Selector);
  IMPL_MEM_OBJ(Type_Selector);
  IMPL_MEM_OBJ(Class_Selector);
  IMPL_MEM_OBJ(Id_Selector);
  IMPL_MEM_OBJ(Attribute_Selector);
  IMPL_MEM_OBJ(Pseudo_Selector);
  IMPL_MEM_OBJ(Wrapped_Selector);
  IMPL_MEM_OBJ(Compound_Selector);
  IMPL_MEM_OBJ(Complex_Selector);
  IMPL_MEM_OBJ(Selector_List);

  // ###########################################################################
  // Implement compare, order and hashing operations for AST Nodes
  // ###########################################################################

  struct HashNodes {
    template <class T>
    size_t operator() (const T& ex) const {
      return ex.isNull() ? 0 : ex->hash();
    }
  };
  template <class T>
  bool OrderFunction(const T& lhs, const T& rhs) {
      return !lhs.isNull() && !rhs.isNull() && *lhs < *rhs;
  };
  struct OrderNodes {
    template <class T>
    bool operator() (const T& lhs, const T& rhs) const {
      return OrderFunction<T>(lhs, rhs);
    }
  };
  template <class T>
  bool CompareFunction(const T& lhs, const T& rhs) {
      // code around sass logic issue. 1px == 1 is true
      // but both items are still different keys in maps
      if (dynamic_cast<Number*>(lhs.ptr()))
        if (dynamic_cast<Number*>(rhs.ptr()))
          return lhs->hash() == rhs->hash();
      return !lhs.isNull() && !rhs.isNull() && *lhs == *rhs;
  }
  struct CompareNodes {
    template <class T>
    bool operator() (const T& lhs, const T& rhs) const {
      return CompareFunction<T>(lhs, rhs);
    }
  };

  struct HashPtr {
    template <class T>
    size_t operator() (const T *ref) const {
      return ref->hash();
    }
  };
  struct ComparePtrs {
    template <class T>
    bool operator() (const T *lhs, const T *rhs) const {
      return *lhs == *rhs;
    }
  };

  // ###########################################################################
  // some often used typedefs
  // ###########################################################################

  typedef std::unordered_map<
    Expression_Obj, // key
    Expression_Obj, // value
    HashNodes, // hasher
    CompareNodes // compare
  > ExpressionMap;
  typedef std::unordered_set<
    Expression_Obj, // value
    HashNodes, // hasher
    CompareNodes // compare
  > ExpressionSet;

  typedef std::string SubSetMapKey;
  typedef std::vector<std::string> SubSetMapKeys;

  typedef std::pair<Complex_Selector_Obj, Compound_Selector_Obj> SubSetMapPair;
  typedef std::pair<Compound_Selector_Obj, Complex_Selector_Obj> SubSetMapLookup;
  typedef std::vector<SubSetMapPair> SubSetMapPairs;
  typedef std::vector<SubSetMapLookup> SubSetMapLookups;

  typedef std::pair<Complex_Selector_Obj, SubSetMapPairs> SubSetMapResult;
  typedef std::vector<SubSetMapResult> SubSetMapResults;

  typedef std::set<Selector_Obj, OrderNodes> SelectorSet;

  typedef std::deque<Complex_Selector_Obj> ComplexSelectorDeque;
  typedef std::set<Simple_Selector_Obj, OrderNodes> SimpleSelectorSet;
  typedef std::set<Complex_Selector_Obj, OrderNodes> ComplexSelectorSet;
  typedef std::set<Compound_Selector_Obj, OrderNodes> CompoundSelectorSet;
  typedef std::unordered_set<Simple_Selector_Obj, HashNodes, CompareNodes> SimpleSelectorDict;

  typedef std::vector<Block_Ptr> BlockStack;
  typedef std::vector<Sass_Callee> CalleeStack;
  typedef std::vector<AST_Node_Obj> CallStack;
  typedef std::vector<Media_Block_Ptr> MediaStack;
  typedef std::vector<Selector_List_Obj> SelectorStack;
  typedef std::vector<Sass_Import_Entry> ImporterStack;

  // only to switch implementations for testing
  #define environment_map std::map

  // ###########################################################################
  // explicit type conversion functions
  // ###########################################################################

  template<class T>
  T* Cast(AST_Node* ptr);

  template<class T>
  const T* Cast(const AST_Node* ptr);

  // sometimes you know the class you want to cast to is final
  // in this case a simple typeid check is faster and safe to use

  #define DECLARE_BASE_CAST(T) \
  template<> T* Cast(AST_Node* ptr); \
  template<> const T* Cast(const AST_Node* ptr); \

  // ###########################################################################
  // implement specialization for final classes
  // ###########################################################################

  DECLARE_BASE_CAST(AST_Node)
  DECLARE_BASE_CAST(Expression)
  DECLARE_BASE_CAST(Statement)
  DECLARE_BASE_CAST(Has_Block)
  DECLARE_BASE_CAST(PreValue)
  DECLARE_BASE_CAST(Value)
  DECLARE_BASE_CAST(List)
  DECLARE_BASE_CAST(Color)
  DECLARE_BASE_CAST(String)
  DECLARE_BASE_CAST(String_Constant)
  DECLARE_BASE_CAST(Supports_Condition)
  DECLARE_BASE_CAST(Selector)
  DECLARE_BASE_CAST(Simple_Selector)

}

#endif

namespace Sass {

  namespace File {

    // return the current directory
    // always with forward slashes
    std::string get_cwd();

    // test if path exists and is a file
    bool file_exists(const std::string& file);

    // return if given path is absolute
    // works with *nix and windows paths
    bool is_absolute_path(const std::string& path);

    // return only the directory part of path
    std::string dir_name(const std::string& path);

    // return only the filename part of path
    std::string base_name(const std::string&);

    // do a locigal clean up of the path
    // no physical check on the filesystem
    std::string make_canonical_path (std::string path);

    // join two path segments cleanly together
    // but only if right side is not absolute yet
    std::string join_paths(std::string root, std::string name);

    // if the relative path is outside of the cwd we want want to
    // show the absolute path in console messages
    std::string path_for_console(const std::string& rel_path, const std::string& abs_path, const std::string& orig_path);

    // create an absolute path by resolving relative paths with cwd
    std::string rel2abs(const std::string& path, const std::string& base = ".", const std::string& cwd = get_cwd());

    // create a path that is relative to the given base directory
    // path and base will first be resolved against cwd to make them absolute
    std::string abs2rel(const std::string& path, const std::string& base = ".", const std::string& cwd = get_cwd());

    // helper function to resolve a filename
    // searching without variations in all paths
    std::string find_file(const std::string& file, struct Sass_Compiler* options);
    std::string find_file(const std::string& file, const std::vector<std::string> paths);

    // helper function to resolve a include filename
    // this has the original resolve logic for sass include
    std::string find_include(const std::string& file, const std::vector<std::string> paths);

    // split a path string delimited by semicolons or colons (OS dependent)
    std::vector<std::string> split_path_list(const char* paths);

    // try to load the given filename
    // returned memory must be freed
    // will auto convert .sass files
    char* read_file(const std::string& file);

  }

  // requested import
  class Importer {
    public:
      // requested import path
      std::string imp_path;
      // parent context path
      std::string ctx_path;
      // base derived from context path
      // this really just acts as a cache
      std::string base_path;
    public:
      Importer(std::string imp_path, std::string ctx_path)
      : imp_path(File::make_canonical_path(imp_path)),
        ctx_path(File::make_canonical_path(ctx_path)),
        base_path(File::dir_name(ctx_path))
      { }
  };

  // a resolved include (final import)
  class Include : public Importer {
    public:
      // resolved absolute path
      std::string abs_path;
    public:
      Include(const Importer& imp, std::string abs_path)
      : Importer(imp), abs_path(abs_path)
      { }
  };

  // a loaded resource
  class Resource {
    public:
      // the file contents
      char* contents;
      // conected sourcemap
      char* srcmap;
    public:
      Resource(char* contents, char* srcmap)
      : contents(contents), srcmap(srcmap)
      { }
  };

  // parsed stylesheet from loaded resource
  class StyleSheet : public Resource {
    public:
      // parsed root block
      Block_Obj root;
    public:
      StyleSheet(const Resource& res, Block_Obj root)
      : Resource(res), root(root)
      { }
  };

  namespace File {

    static std::vector<std::string> defaultExtensions = { ".scss", ".sass", ".css" };

    std::vector<Include> resolve_includes(const std::string& root, const std::string& file,
      const std::vector<std::string>& exts = defaultExtensions);

  }

}

#endif

namespace Sass {

  struct Backtrace {

    ParserState pstate;
    std::string caller;

    Backtrace(ParserState pstate, std::string c = "")
    : pstate(pstate),
      caller(c)
    { }

  };

  typedef std::vector<Backtrace> Backtraces;

  const std::string traces_to_string(Backtraces traces, std::string indent = "\t");

}

#endif

namespace Sass {

  Value_Ptr c2ast(union Sass_Value* v, Backtraces traces, ParserState pstate);

}

#endif
#ifndef SASS_EXTEND_H
#define SASS_EXTEND_H


#ifndef SASS_AST_H
#define SASS_AST_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.
// must be the first include in all compile units
#ifndef SASS_SASS_H
#define SASS_SASS_H

// undefine extensions macro to tell sys includes
// that we do not want any macros to be exported
// mainly fixes an issue on SmartOS (SEC macro)
#undef __EXTENSIONS__

#ifdef _MSC_VER
#pragma warning(disable : 4005)
#endif

// aplies to MSVC and MinGW
#ifdef _WIN32
// we do not want the ERROR macro
# ifndef NOGDI
#  define NOGDI
# endif
// we do not want the min/max macro
# ifndef NOMINMAX
#  define NOMINMAX
# endif
// we do not want the IN/OUT macro
# ifndef _NO_W32_PSEUDO_MODIFIERS
#  define _NO_W32_PSEUDO_MODIFIERS
# endif
#endif


// should we be case insensitive
// when dealing with files or paths
#ifndef FS_CASE_SENSITIVE
# ifdef _WIN32
#  define FS_CASE_SENSITIVE 0
# else
#  define FS_CASE_SENSITIVE 1
# endif
#endif

// path separation char
#ifndef PATH_SEP
# ifdef _WIN32
#  define PATH_SEP ';'
# else
#  define PATH_SEP ':'
# endif
#endif


// include C-API header

// For C++ helper

// output behaviours
namespace Sass {

  // create some C++ aliases for the most used options
  const static Sass_Output_Style NESTED = SASS_STYLE_NESTED;
  const static Sass_Output_Style COMPACT = SASS_STYLE_COMPACT;
  const static Sass_Output_Style EXPANDED = SASS_STYLE_EXPANDED;
  const static Sass_Output_Style COMPRESSED = SASS_STYLE_COMPRESSED;
  // only used internal to trigger ruby inspect behavior
  const static Sass_Output_Style INSPECT = SASS_STYLE_INSPECT;
  const static Sass_Output_Style TO_SASS = SASS_STYLE_TO_SASS;

  // helper to aid dreaded MSVC debug mode
  // see implementation for more details
  char* sass_copy_string(std::string str);

}

// input behaviours
enum Sass_Input_Style {
  SASS_CONTEXT_NULL,
  SASS_CONTEXT_FILE,
  SASS_CONTEXT_DATA,
  SASS_CONTEXT_FOLDER
};

// simple linked list
struct string_list {
  string_list* next;
  char* string;
};

// sass config options structure
struct Sass_Inspect_Options {

  // Output style for the generated css code
  // A value from above SASS_STYLE_* constants
  enum Sass_Output_Style output_style;

  // Precision for fractional numbers
  int precision;

  // initialization list (constructor with defaults)
  Sass_Inspect_Options(Sass_Output_Style style = Sass::NESTED,
                       int precision = 10)
  : output_style(style), precision(precision)
  { }

};

// sass config options structure
struct Sass_Output_Options : Sass_Inspect_Options {

  // String to be used for indentation
  const char* indent;
  // String to be used to for line feeds
  const char* linefeed;

  // Emit comments in the generated CSS indicating
  // the corresponding source line.
  bool source_comments;

  // initialization list (constructor with defaults)
  Sass_Output_Options(struct Sass_Inspect_Options opt,
                      const char* indent = "  ",
                      const char* linefeed = "\n",
                      bool source_comments = false)
  : Sass_Inspect_Options(opt),
    indent(indent), linefeed(linefeed),
    source_comments(source_comments)
  { }

  // initialization list (constructor with defaults)
  Sass_Output_Options(Sass_Output_Style style = Sass::NESTED,
                      int precision = 10,
                      const char* indent = "  ",
                      const char* linefeed = "\n",
                      bool source_comments = false)
  : Sass_Inspect_Options(style, precision),
    indent(indent), linefeed(linefeed),
    source_comments(source_comments)
  { }

};

#endif


#ifndef SASS_UTIL_H
#define SASS_UTIL_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.


#include <assert.h>
#include <math.h>

#define SASS_ASSERT(cond, msg) assert(cond && msg)

namespace Sass {

  template <typename T>
  T clip(const T& n, const T& lower, const T& upper) {
    return std::max(lower, std::min(n, upper));
  }

  template <typename T>
  T absmod(const T& n, const T& r) {
    T m = std::fmod(n, r);
    if (m < 0.0) m += r;
    return m;
  }

  double round(double val, size_t precision = 0);
  double sass_strtod(const char* str);
  const char* safe_str(const char *, const char* = "");
  void free_string_array(char **);
  char **copy_strings(const std::vector<std::string>&, char ***, int = 0);
  std::string read_css_string(const std::string& str, bool css = true);
  std::string evacuate_escapes(const std::string& str);
  std::string string_to_output(const std::string& str);
  std::string comment_to_string(const std::string& text);
  std::string read_hex_escapes(const std::string& str);
  std::string escape_string(const std::string& str);
  void newline_to_space(std::string& str);

  std::string quote(const std::string&, char q = 0);
  std::string unquote(const std::string&, char* q = 0, bool keep_utf8_sequences = false, bool strict = true);
  char detect_best_quotemark(const char* s, char qm = '"');

  bool is_hex_doublet(double n);
  bool is_color_doublet(double r, double g, double b);

  bool peek_linefeed(const char* start);

  // Returns true iff `elements` ⊆ `container`.
  template <typename C, typename T>
  bool contains_all(C container, T elements) {
    for (const auto &el : elements) {
      if (container.find(el) == container.end()) return false;
    }
    return true;
  }

  // C++20 `starts_with` equivalent.
  // See https://en.cppreference.com/w/cpp/string/basic_string/starts_with
  inline bool starts_with(const std::string& str, const char* prefix, size_t prefix_len) {
    return str.compare(0, prefix_len, prefix) == 0;
  }

  inline bool starts_with(const std::string& str, const char* prefix) {
    return starts_with(str, prefix, std::strlen(prefix));
  }

  // C++20 `ends_with` equivalent.
  // See https://en.cppreference.com/w/cpp/string/basic_string/ends_with
  inline bool ends_with(const std::string& str, const std::string& suffix) {
    return suffix.size() <= str.size() && std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
  }

  inline bool ends_with(const std::string& str, const char* suffix, size_t suffix_len) {
    if (suffix_len > str.size()) return false;
    const char* suffix_it = suffix + suffix_len;
    const char* str_it = str.c_str() + str.size();
    while (suffix_it != suffix) if (*(--suffix_it) != *(--str_it)) return false;
    return true;
  }

  inline bool ends_with(const std::string& str, const char* suffix) {
    return ends_with(str, suffix, std::strlen(suffix));
  }

  namespace Util {

    std::string rtrim(const std::string& str);

    std::string normalize_underscores(const std::string& str);
    std::string normalize_decimals(const std::string& str);

    bool isPrintable(Ruleset_Ptr r, Sass_Output_Style style = NESTED);
    bool isPrintable(Supports_Block_Ptr r, Sass_Output_Style style = NESTED);
    bool isPrintable(Media_Block_Ptr r, Sass_Output_Style style = NESTED);
    bool isPrintable(Comment_Ptr b, Sass_Output_Style style = NESTED);
    bool isPrintable(Block_Obj b, Sass_Output_Style style = NESTED);
    bool isPrintable(String_Constant_Ptr s, Sass_Output_Style style = NESTED);
    bool isPrintable(String_Quoted_Ptr s, Sass_Output_Style style = NESTED);
    bool isPrintable(Declaration_Ptr d, Sass_Output_Style style = NESTED);
    bool isAscii(const char chr);

  }
}
#endif
#ifndef SASS_UNITS_H
#define SASS_UNITS_H

#include <cmath>

namespace Sass {

  const double PI = std::acos(-1);

  enum UnitClass {
    LENGTH = 0x000,
    ANGLE = 0x100,
    TIME = 0x200,
    FREQUENCY = 0x300,
    RESOLUTION = 0x400,
    INCOMMENSURABLE = 0x500
  };

  enum UnitType {

    // size units
    IN = UnitClass::LENGTH,
    CM,
    PC,
    MM,
    PT,
    PX,

    // angle units
    DEG = ANGLE,
    GRAD,
    RAD,
    TURN,

    // time units
    SEC = TIME,
    MSEC,

    // frequency units
    HERTZ = FREQUENCY,
    KHERTZ,

    // resolutions units
    DPI = RESOLUTION,
    DPCM,
    DPPX,

    // for unknown units
    UNKNOWN = INCOMMENSURABLE

  };

  class Units {
  public:
    std::vector<std::string> numerators;
    std::vector<std::string> denominators;
  public:
    // default constructor
    Units() :
      numerators(),
      denominators()
    { }
    // copy constructor
    Units(const Units* ptr) :
      numerators(ptr->numerators),
      denominators(ptr->denominators)
    { }
    // convert to string
    std::string unit() const;
    // get if units are empty
    bool is_unitless() const;
    // return if valid for css
    bool is_valid_css_unit() const;
    // reduce units for output
    // returns conversion factor
    double reduce();
    // normalize units for compare
    // returns conversion factor
    double normalize();
    // compare operations
    bool operator< (const Units& rhs) const;
    bool operator== (const Units& rhs) const;
    bool operator!= (const Units& rhs) const;
    // factor to convert into given units
    double convert_factor(const Units&) const;
  };

  extern const double size_conversion_factors[6][6];
  extern const double angle_conversion_factors[4][4];
  extern const double time_conversion_factors[2][2];
  extern const double frequency_conversion_factors[2][2];
  extern const double resolution_conversion_factors[3][3];

  UnitType get_main_unit(const UnitClass unit);
  enum Sass::UnitType string_to_unit(const std::string&);
  const char* unit_to_string(Sass::UnitType unit);
  enum Sass::UnitClass get_unit_type(Sass::UnitType unit);
  std::string get_unit_class(Sass::UnitType unit);
  std::string unit_to_class(const std::string&);
  // throws incompatibleUnits exceptions
  double conversion_factor(const std::string&, const std::string&);
  double conversion_factor(UnitType, UnitType, UnitClass, UnitClass);
  double convert_units(const std::string&, const std::string&, int&, int&);

}

#endif
#ifndef SASS_CONTEXT_H
#define SASS_CONTEXT_H


#define BUFFERSIZE 255
// :mode=c++:
/*
encode.h - c++ wrapper for a base64 encoding algorithm

This is part of the libb64 project, and has been placed in the public domain.
For details, see http://sourceforge.net/projects/libb64
*/
#ifndef BASE64_ENCODE_H
#define BASE64_ENCODE_H


namespace base64
{
	extern "C"
	{
/*
cencode.h - c header for a base64 encoding algorithm

This is part of the libb64 project, and has been placed in the public domain.
For details, see http://sourceforge.net/projects/libb64
*/

#ifndef BASE64_CENCODE_H
#define BASE64_CENCODE_H

typedef enum
{
  step_A, step_B, step_C
} base64_encodestep;

typedef struct
{
	base64_encodestep step;
	char result;
	int stepcount;
} base64_encodestate;

void base64_init_encodestate(base64_encodestate* state_in);

char base64_encode_value(char value_in);

int base64_encode_block(const char* plaintext_in, int length_in, char* code_out, base64_encodestate* state_in);

int base64_encode_blockend(char* code_out, base64_encodestate* state_in);

#endif /* BASE64_CENCODE_H */

	}

	struct encoder
	{
		base64_encodestate _state;
		int _buffersize;

		encoder(int buffersize_in = BUFFERSIZE)
		: _buffersize(buffersize_in)
		{
			base64_init_encodestate(&_state);
		}

		int encode(char value_in)
		{
			return base64_encode_value(value_in);
		}

		int encode(const char* code_in, const int length_in, char* plaintext_out)
		{
			return base64_encode_block(code_in, length_in, plaintext_out, &_state);
		}

		int encode_end(char* plaintext_out)
		{
			return base64_encode_blockend(plaintext_out, &_state);
		}

		void encode(std::istream& istream_in, std::ostream& ostream_in)
		{
			base64_init_encodestate(&_state);
			//
			const int N = _buffersize;
			char* plaintext = new char[N];
			char* code = new char[2*N];
			int plainlength;
			int codelength;

			do
			{
				istream_in.read(plaintext, N);
				plainlength = static_cast<int>(istream_in.gcount());
				//
				codelength = encode(plaintext, plainlength, code);
				ostream_in.write(code, codelength);
			}
			while (istream_in.good() && plainlength > 0);

			codelength = encode_end(code);
			ostream_in.write(code, codelength);
			//
			base64_init_encodestate(&_state);

			delete [] code;
			delete [] plaintext;
		}
	};

} // namespace base64

#endif // BASE64_ENCODE_H


#ifndef SASS_KWD_ARG_MACROS_H
#define SASS_KWD_ARG_MACROS_H

// Example usage:
// KWD_ARG_SET(Args) {
//   KWD_ARG(Args, string, foo);
//   KWD_ARG(Args, int, bar);
//   ...
// };
//
// ... and later ...
//
// something(Args().foo("hey").bar(3));

#define KWD_ARG_SET(set_name) class set_name

#define KWD_ARG(set_name, type, name) \
private: \
  type name##_; \
public: \
  set_name& name(type name##__) { \
    name##_ = name##__; \
    return *this; \
  } \
  type name() { return name##_; } \
private:

#endif
#ifndef SASS_SASS_CONTEXT_H
#define SASS_SASS_CONTEXT_H


// sass config options structure
struct Sass_Options : Sass_Output_Options {

  // embed sourceMappingUrl as data uri
  bool source_map_embed;

  // embed include contents in maps
  bool source_map_contents;

  // create file urls for sources
  bool source_map_file_urls;

  // Disable sourceMappingUrl in css output
  bool omit_source_map_url;

  // Treat source_string as sass (as opposed to scss)
  bool is_indented_syntax_src;

  // The input path is used for source map
  // generation. It can be used to define
  // something with string compilation or to
  // overload the input file path. It is
  // set to "stdin" for data contexts and
  // to the input file on file contexts.
  char* input_path;

  // The output path is used for source map
  // generation. LibSass will not write to
  // this file, it is just used to create
  // information in source-maps etc.
  char* output_path;

  // Colon-separated list of paths
  // Semicolon-separated on Windows
  // Maybe use array interface instead?
  char* include_path;
  char* plugin_path;

  // Include paths (linked string list)
  struct string_list* include_paths;
  // Plugin paths (linked string list)
  struct string_list* plugin_paths;

  // Path to source map file
  // Enables source map generation
  // Used to create sourceMappingUrl
  char* source_map_file;

  // Directly inserted in source maps
  char* source_map_root;

  // Custom functions that can be called from sccs code
  Sass_Function_List c_functions;

  // List of custom importers
  Sass_Importer_List c_importers;

  // List of custom headers
  Sass_Importer_List c_headers;

};


// base for all contexts
struct Sass_Context : Sass_Options
{

  // store context type info
  enum Sass_Input_Style type;

  // generated output data
  char* output_string;

  // generated source map json
  char* source_map_string;

  // error status
  int error_status;
  char* error_json;
  char* error_text;
  char* error_message;
  // error position
  char* error_file;
  size_t error_line;
  size_t error_column;
  const char* error_src;

  // report imported files
  char** included_files;

};

// struct for file compilation
struct Sass_File_Context : Sass_Context {

  // no additional fields required
  // input_path is already on options

};

// struct for data compilation
struct Sass_Data_Context : Sass_Context {

  // provided source string
  char* source_string;
  char* srcmap_string;

};

// link c and cpp context
struct Sass_Compiler {
  // progress status
  Sass_Compiler_State state;
  // original c context
  Sass_Context* c_ctx;
  // Sass::Context
  Sass::Context* cpp_ctx;
  // Sass::Block
  Sass::Block_Obj root;
};

#endif
#ifndef SASS_ENVIRONMENT_H
#define SASS_ENVIRONMENT_H

#ifndef SASS_AST_DEF_MACROS_H
#define SASS_AST_DEF_MACROS_H

// Helper class to switch a flag and revert once we go out of scope
template <class T>
class LocalOption {
  private:
    T* var; // pointer to original variable
    T orig; // copy of the original option
  public:
    LocalOption(T& var)
    {
      this->var = &var;
      this->orig = var;
    }
    LocalOption(T& var, T orig)
    {
      this->var = &var;
      this->orig = var;
      *(this->var) = orig;
    }
    void reset()
    {
      *(this->var) = this->orig;
    }
    ~LocalOption() {
      *(this->var) = this->orig;
    }
};

#define LOCAL_FLAG(name,opt) LocalOption<bool> flag_##name(name, opt)
#define LOCAL_COUNT(name,opt) LocalOption<size_t> cnt_##name(name, opt)

#define NESTING_GUARD(name) \
  LocalOption<size_t> cnt_##name(name, name + 1); \
  if (name > MAX_NESTING) throw Exception::NestingLimitError(pstate, traces); \

#define ADD_PROPERTY(type, name)\
protected:\
  type name##_;\
public:\
  type name() const        { return name##_; }\
  type name(type name##__) { return name##_ = name##__; }\
private:

#define HASH_PROPERTY(type, name)\
protected:\
  type name##_;\
public:\
  type name() const        { return name##_; }\
  type name(type name##__) { hash_ = 0; return name##_ = name##__; }\
private:

#define ADD_CONSTREF(type, name) \
protected: \
  type name##_; \
public: \
  const type& name() const { return name##_; } \
  void name(type name##__) { name##_ = name##__; } \
private:

#define HASH_CONSTREF(type, name) \
protected: \
  type name##_; \
public: \
  const type& name() const { return name##_; } \
  void name(type name##__) { hash_ = 0; name##_ = name##__; } \
private:

#ifdef DEBUG_SHARED_PTR

#define ATTACH_ABSTRACT_AST_OPERATIONS(klass) \
  virtual klass##_Ptr copy(std::string, size_t) const = 0; \
  virtual klass##_Ptr clone(std::string, size_t) const = 0; \

#define ATTACH_VIRTUAL_AST_OPERATIONS(klass) \
  klass(const klass* ptr); \
  virtual klass##_Ptr copy(std::string, size_t) const override = 0; \
  virtual klass##_Ptr clone(std::string, size_t) const override = 0; \

#define ATTACH_AST_OPERATIONS(klass) \
  klass(const klass* ptr); \
  virtual klass##_Ptr copy(std::string, size_t) const override; \
  virtual klass##_Ptr clone(std::string, size_t) const override; \

#else

#define ATTACH_ABSTRACT_AST_OPERATIONS(klass) \
  virtual klass##_Ptr copy() const = 0; \
  virtual klass##_Ptr clone() const = 0; \

#define ATTACH_VIRTUAL_AST_OPERATIONS(klass) \
  klass(const klass* ptr); \
  virtual klass##_Ptr copy() const override = 0; \
  virtual klass##_Ptr clone() const override = 0; \

#define ATTACH_AST_OPERATIONS(klass) \
  klass(const klass* ptr); \
  virtual klass##_Ptr copy() const override; \
  virtual klass##_Ptr clone() const override; \

#endif

#ifdef DEBUG_SHARED_PTR

  #define IMPLEMENT_AST_OPERATORS(klass) \
    klass##_Ptr klass::copy(std::string file, size_t line) const { \
      klass##_Ptr cpy = new klass(this); \
      cpy->trace(file, line); \
      return cpy; \
    } \
    klass##_Ptr klass::clone(std::string file, size_t line) const { \
      klass##_Ptr cpy = copy(file, line); \
      cpy->cloneChildren(); \
      return cpy; \
    } \

#else

  #define IMPLEMENT_AST_OPERATORS(klass) \
    klass##_Ptr klass::copy() const { \
      return new klass(this); \
    } \
    klass##_Ptr klass::clone() const { \
      klass##_Ptr cpy = copy(); \
      cpy->cloneChildren(); \
      return cpy; \
    } \

#endif

#endif

namespace Sass {

  // this defeats the whole purpose of environment being templatable!!
  typedef environment_map<std::string, AST_Node_Obj>::iterator EnvIter;

  class EnvResult {
    public:
      EnvIter it;
      bool found;
    public:
      EnvResult(EnvIter it, bool found)
      : it(it), found(found) {}
  };

  template <typename T>
  class Environment {
    // TODO: test with map
    environment_map<std::string, T> local_frame_;
    ADD_PROPERTY(Environment*, parent)
    ADD_PROPERTY(bool, is_shadow)

  public:
    Environment(bool is_shadow = false);
    Environment(Environment* env, bool is_shadow = false);
    Environment(Environment& env, bool is_shadow = false);

    // link parent to create a stack
    void link(Environment& env);
    void link(Environment* env);

    // this is used to find the global frame
    // which is the second last on the stack
    bool is_lexical() const;

    // only match the real root scope
    // there is still a parent around
    // not sure what it is actually use for
    // I guess we store functions etc. there
    bool is_global() const;

    // scope operates on the current frame

    environment_map<std::string, T>& local_frame();

    bool has_local(const std::string& key) const;

    EnvResult find_local(const std::string& key);

    T& get_local(const std::string& key);

    // set variable on the current frame
    void set_local(const std::string& key, const T& val);
    void set_local(const std::string& key, T&& val);

    void del_local(const std::string& key);

    // global operates on the global frame
    // which is the second last on the stack
    Environment* global_env();
    // get the env where the variable already exists
    // if it does not yet exist, we return current env
    Environment* lexical_env(const std::string& key);

    bool has_global(const std::string& key);

    T& get_global(const std::string& key);

    // set a variable on the global frame
    void set_global(const std::string& key, const T& val);
    void set_global(const std::string& key, T&& val);

    void del_global(const std::string& key);

    // see if we have a lexical variable
    // move down the stack but stop before we
    // reach the global frame (is not included)
    bool has_lexical(const std::string& key) const;

    // see if we have a lexical we could update
    // either update already existing lexical value
    // or we create a new one on the current frame
    void set_lexical(const std::string& key, T&& val);
    void set_lexical(const std::string& key, const T& val);

    // look on the full stack for key
    // include all scopes available
    bool has(const std::string& key) const;

    // look on the full stack for key
    // include all scopes available
    T& get(const std::string& key);

    // look on the full stack for key
    // include all scopes available
    EnvResult find(const std::string& key);

    // use array access for getter and setter functions
    T& operator[](const std::string& key);

    #ifdef DEBUG
    size_t print(std::string prefix = "");
    #endif

  };

  // define typedef for our use case
  typedef Environment<AST_Node_Obj> Env;
  typedef std::vector<Env*> EnvStack;

}

#endif
#ifndef SASS_SOURCE_MAP_H
#define SASS_SOURCE_MAP_H


#ifndef SASS_BASE64VLQ_H
#define SASS_BASE64VLQ_H


namespace Sass {

  class Base64VLQ {

  public:

    std::string encode(const int number) const;

  private:

    char base64_encode(const int number) const;

    int to_vlq_signed(const int number) const;

    static const char* CHARACTERS;

    static const int VLQ_BASE_SHIFT;
    static const int VLQ_BASE;
    static const int VLQ_BASE_MASK;
    static const int VLQ_CONTINUATION_BIT;
  };

}

#endif

#define VECTOR_PUSH(vec, ins) vec.insert(vec.end(), ins.begin(), ins.end())
#define VECTOR_UNSHIFT(vec, ins) vec.insert(vec.begin(), ins.begin(), ins.end())

namespace Sass {

  class Context;
  class OutputBuffer;

  class SourceMap {

  public:
    std::vector<size_t> source_index;
    SourceMap();
    SourceMap(const std::string& file);

    void append(const Offset& offset);
    void prepend(const Offset& offset);
    void append(const OutputBuffer& out);
    void prepend(const OutputBuffer& out);
    void add_open_mapping(AST_Node_Ptr_Const node);
    void add_close_mapping(AST_Node_Ptr_Const node);

    std::string render_srcmap(Context &ctx);
    ParserState remap(const ParserState& pstate);

  private:

    std::string serialize_mappings();

    std::vector<Mapping> mappings;
    Position current_position;
public:
    std::string file;
private:
    Base64VLQ base64vlq;
  };

  class OutputBuffer {
    public:
      OutputBuffer(void)
      : buffer(""),
        smap()
      { }
    public:
      std::string buffer;
      SourceMap smap;
  };

}

#endif
#ifndef SASS_SUBSET_MAP_H
#define SASS_SUBSET_MAP_H

#include <iterator>



// #include <iostream>
// #include <sstream>
// template<typename T>
// std::string vector_to_string(std::vector<T> v)
// {
//   std::stringstream buffer;
//   buffer << "[";

//   if (!v.empty())
//   {  buffer << v[0]; }
//   else
//   { buffer << "]"; }

//   if (v.size() == 1)
//   { buffer << "]"; }
//   else
//   {
//     for (size_t i = 1, S = v.size(); i < S; ++i) buffer << ", " << v[i];
//     buffer << "]";
//   }

//   return buffer.str();
// }

// template<typename T>
// std::string set_to_string(set<T> v)
// {
//   std::stringstream buffer;
//   buffer << "[";
//   typename std::set<T>::iterator i = v.begin();
//   if (!v.empty())
//   {  buffer << *i; }
//   else
//   { buffer << "]"; }

//   if (v.size() == 1)
//   { buffer << "]"; }
//   else
//   {
//     for (++i; i != v.end(); ++i) buffer << ", " << *i;
//     buffer << "]";
//   }

//   return buffer.str();
// }

namespace Sass {

  class Subset_Map {
  private:
    std::vector<SubSetMapPair> values_;
    std::map<Simple_Selector_Obj, std::vector<std::pair<Compound_Selector_Obj, size_t> >, OrderNodes > hash_;
  public:
    void put(const Compound_Selector_Obj& sel, const SubSetMapPair& value);
    std::vector<SubSetMapPair> get_kv(const Compound_Selector_Obj& s);
    std::vector<SubSetMapPair> get_v(const Compound_Selector_Obj& s);
    bool empty() { return values_.empty(); }
    void clear() { values_.clear(); hash_.clear(); }
    const std::vector<SubSetMapPair> values(void) { return values_; }
  };

}

#endif
#ifndef SASS_OUTPUT_H
#define SASS_OUTPUT_H


#ifndef SASS_INSPECT_H
#define SASS_INSPECT_H

#ifndef SASS_OPERATION_H
#define SASS_OPERATION_H

// base classes to implement curiously recurring template pattern (CRTP)
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern

#include <stdexcept>


namespace Sass {

  #define ATTACH_ABSTRACT_CRTP_PERFORM_METHODS()\
    virtual void perform(Operation<void>* op) = 0; \
    virtual Value_Ptr perform(Operation<Value_Ptr>* op) = 0; \
    virtual std::string perform(Operation<std::string>* op) = 0; \
    virtual AST_Node_Ptr perform(Operation<AST_Node_Ptr>* op) = 0; \
    virtual Selector_Ptr perform(Operation<Selector_Ptr>* op) = 0; \
    virtual Statement_Ptr perform(Operation<Statement_Ptr>* op) = 0; \
    virtual Expression_Ptr perform(Operation<Expression_Ptr>* op) = 0; \
    virtual union Sass_Value* perform(Operation<union Sass_Value*>* op) = 0; \
    virtual Supports_Condition_Ptr perform(Operation<Supports_Condition_Ptr>* op) = 0; \

  // you must add operators to every class
  // ensures `this` of actual instance type
  // we therefore call the specific operator
  // they are virtual so most specific is used
  #define ATTACH_CRTP_PERFORM_METHODS()\
    virtual void perform(Operation<void>* op) override { return (*op)(this); } \
    virtual Value_Ptr perform(Operation<Value_Ptr>* op) override { return (*op)(this); } \
    virtual std::string perform(Operation<std::string>* op) override { return (*op)(this); } \
    virtual AST_Node_Ptr perform(Operation<AST_Node_Ptr>* op) override { return (*op)(this); } \
    virtual Selector_Ptr perform(Operation<Selector_Ptr>* op) override { return (*op)(this); } \
    virtual Statement_Ptr perform(Operation<Statement_Ptr>* op) override { return (*op)(this); } \
    virtual Expression_Ptr perform(Operation<Expression_Ptr>* op) override { return (*op)(this); } \
    virtual union Sass_Value* perform(Operation<union Sass_Value*>* op) override { return (*op)(this); } \
    virtual Supports_Condition_Ptr perform(Operation<Supports_Condition_Ptr>* op) override { return (*op)(this); } \

  template<typename T>
  class Operation {
  public:
    virtual T operator()(AST_Node_Ptr x)               = 0;
    // statements
    virtual T operator()(Block_Ptr x)                  = 0;
    virtual T operator()(Ruleset_Ptr x)                = 0;
    virtual T operator()(Bubble_Ptr x)                 = 0;
    virtual T operator()(Trace_Ptr x)                  = 0;
    virtual T operator()(Supports_Block_Ptr x)         = 0;
    virtual T operator()(Media_Block_Ptr x)            = 0;
    virtual T operator()(At_Root_Block_Ptr x)          = 0;
    virtual T operator()(Directive_Ptr x)              = 0;
    virtual T operator()(Keyframe_Rule_Ptr x)          = 0;
    virtual T operator()(Declaration_Ptr x)            = 0;
    virtual T operator()(Assignment_Ptr x)             = 0;
    virtual T operator()(Import_Ptr x)                 = 0;
    virtual T operator()(Import_Stub_Ptr x)            = 0;
    virtual T operator()(Warning_Ptr x)                = 0;
    virtual T operator()(Error_Ptr x)                  = 0;
    virtual T operator()(Debug_Ptr x)                  = 0;
    virtual T operator()(Comment_Ptr x)                = 0;
    virtual T operator()(If_Ptr x)                     = 0;
    virtual T operator()(For_Ptr x)                    = 0;
    virtual T operator()(Each_Ptr x)                   = 0;
    virtual T operator()(While_Ptr x)                  = 0;
    virtual T operator()(Return_Ptr x)                 = 0;
    virtual T operator()(Content_Ptr x)                = 0;
    virtual T operator()(Extension_Ptr x)              = 0;
    virtual T operator()(Definition_Ptr x)             = 0;
    virtual T operator()(Mixin_Call_Ptr x)             = 0;
    // expressions
    virtual T operator()(Null_Ptr x)                   = 0;
    virtual T operator()(List_Ptr x)                   = 0;
    virtual T operator()(Map_Ptr x)                    = 0;
    virtual T operator()(Function_Ptr x)               = 0;
    virtual T operator()(Binary_Expression_Ptr x)      = 0;
    virtual T operator()(Unary_Expression_Ptr x)       = 0;
    virtual T operator()(Function_Call_Ptr x)          = 0;
    virtual T operator()(Custom_Warning_Ptr x)         = 0;
    virtual T operator()(Custom_Error_Ptr x)           = 0;
    virtual T operator()(Variable_Ptr x)               = 0;
    virtual T operator()(Number_Ptr x)                 = 0;
    virtual T operator()(Color_Ptr x)                  = 0;
    virtual T operator()(Color_RGBA_Ptr x)             = 0;
    virtual T operator()(Color_HSLA_Ptr x)             = 0;
    virtual T operator()(Boolean_Ptr x)                = 0;
    virtual T operator()(String_Schema_Ptr x)          = 0;
    virtual T operator()(String_Quoted_Ptr x)          = 0;
    virtual T operator()(String_Constant_Ptr x)        = 0;
    virtual T operator()(Supports_Condition_Ptr x)     = 0;
    virtual T operator()(Supports_Operator_Ptr x)      = 0;
    virtual T operator()(Supports_Negation_Ptr x)      = 0;
    virtual T operator()(Supports_Declaration_Ptr x)   = 0;
    virtual T operator()(Supports_Interpolation_Ptr x) = 0;
    virtual T operator()(Media_Query_Ptr x)            = 0;
    virtual T operator()(Media_Query_Expression_Ptr x) = 0;
    virtual T operator()(At_Root_Query_Ptr x)          = 0;
    virtual T operator()(Parent_Selector_Ptr x)        = 0;
    virtual T operator()(Parent_Reference_Ptr x)        = 0;
    // parameters and arguments
    virtual T operator()(Parameter_Ptr x)              = 0;
    virtual T operator()(Parameters_Ptr x)             = 0;
    virtual T operator()(Argument_Ptr x)               = 0;
    virtual T operator()(Arguments_Ptr x)              = 0;
    // selectors
    virtual T operator()(Selector_Schema_Ptr x)        = 0;
    virtual T operator()(Placeholder_Selector_Ptr x)   = 0;
    virtual T operator()(Type_Selector_Ptr x)       = 0;
    virtual T operator()(Class_Selector_Ptr x)         = 0;
    virtual T operator()(Id_Selector_Ptr x)            = 0;
    virtual T operator()(Attribute_Selector_Ptr x)     = 0;
    virtual T operator()(Pseudo_Selector_Ptr x)        = 0;
    virtual T operator()(Wrapped_Selector_Ptr x)       = 0;
    virtual T operator()(Compound_Selector_Ptr x)= 0;
    virtual T operator()(Complex_Selector_Ptr x)      = 0;
    virtual T operator()(Selector_List_Ptr x) = 0;
  };

  // example: Operation_CRTP<Expression_Ptr, Eval>
  // T is the base return type of all visitors
  // D is the class derived visitor class
  // normaly you want to implement all operators
  template <typename T, typename D>
  class Operation_CRTP : public Operation<T> {
  public:
    T operator()(AST_Node_Ptr x)               { return static_cast<D*>(this)->fallback(x); }
    // statements
    T operator()(Block_Ptr x)                  { return static_cast<D*>(this)->fallback(x); }
    T operator()(Ruleset_Ptr x)                { return static_cast<D*>(this)->fallback(x); }
    T operator()(Bubble_Ptr x)                 { return static_cast<D*>(this)->fallback(x); }
    T operator()(Trace_Ptr x)                  { return static_cast<D*>(this)->fallback(x); }
    T operator()(Supports_Block_Ptr x)         { return static_cast<D*>(this)->fallback(x); }
    T operator()(Media_Block_Ptr x)            { return static_cast<D*>(this)->fallback(x); }
    T operator()(At_Root_Block_Ptr x)          { return static_cast<D*>(this)->fallback(x); }
    T operator()(Directive_Ptr x)              { return static_cast<D*>(this)->fallback(x); }
    T operator()(Keyframe_Rule_Ptr x)          { return static_cast<D*>(this)->fallback(x); }
    T operator()(Declaration_Ptr x)            { return static_cast<D*>(this)->fallback(x); }
    T operator()(Assignment_Ptr x)             { return static_cast<D*>(this)->fallback(x); }
    T operator()(Import_Ptr x)                 { return static_cast<D*>(this)->fallback(x); }
    T operator()(Import_Stub_Ptr x)            { return static_cast<D*>(this)->fallback(x); }
    T operator()(Warning_Ptr x)                { return static_cast<D*>(this)->fallback(x); }
    T operator()(Error_Ptr x)                  { return static_cast<D*>(this)->fallback(x); }
    T operator()(Debug_Ptr x)                  { return static_cast<D*>(this)->fallback(x); }
    T operator()(Comment_Ptr x)                { return static_cast<D*>(this)->fallback(x); }
    T operator()(If_Ptr x)                     { return static_cast<D*>(this)->fallback(x); }
    T operator()(For_Ptr x)                    { return static_cast<D*>(this)->fallback(x); }
    T operator()(Each_Ptr x)                   { return static_cast<D*>(this)->fallback(x); }
    T operator()(While_Ptr x)                  { return static_cast<D*>(this)->fallback(x); }
    T operator()(Return_Ptr x)                 { return static_cast<D*>(this)->fallback(x); }
    T operator()(Content_Ptr x)                { return static_cast<D*>(this)->fallback(x); }
    T operator()(Extension_Ptr x)              { return static_cast<D*>(this)->fallback(x); }
    T operator()(Definition_Ptr x)             { return static_cast<D*>(this)->fallback(x); }
    T operator()(Mixin_Call_Ptr x)             { return static_cast<D*>(this)->fallback(x); }
    // expressions
    T operator()(Null_Ptr x)                   { return static_cast<D*>(this)->fallback(x); }
    T operator()(List_Ptr x)                   { return static_cast<D*>(this)->fallback(x); }
    T operator()(Map_Ptr x)                    { return static_cast<D*>(this)->fallback(x); }
    T operator()(Function_Ptr x)               { return static_cast<D*>(this)->fallback(x); }
    T operator()(Binary_Expression_Ptr x)      { return static_cast<D*>(this)->fallback(x); }
    T operator()(Unary_Expression_Ptr x)       { return static_cast<D*>(this)->fallback(x); }
    T operator()(Function_Call_Ptr x)          { return static_cast<D*>(this)->fallback(x); }
    T operator()(Custom_Warning_Ptr x)         { return static_cast<D*>(this)->fallback(x); }
    T operator()(Custom_Error_Ptr x)           { return static_cast<D*>(this)->fallback(x); }
    T operator()(Variable_Ptr x)               { return static_cast<D*>(this)->fallback(x); }
    T operator()(Number_Ptr x)                 { return static_cast<D*>(this)->fallback(x); }
    T operator()(Color_Ptr x)                  { return static_cast<D*>(this)->fallback(x); }
    T operator()(Color_RGBA_Ptr x)             { return static_cast<D*>(this)->fallback(x); }
    T operator()(Color_HSLA_Ptr x)             { return static_cast<D*>(this)->fallback(x); }
    T operator()(Boolean_Ptr x)                { return static_cast<D*>(this)->fallback(x); }
    T operator()(String_Schema_Ptr x)          { return static_cast<D*>(this)->fallback(x); }
    T operator()(String_Constant_Ptr x)        { return static_cast<D*>(this)->fallback(x); }
    T operator()(String_Quoted_Ptr x)          { return static_cast<D*>(this)->fallback(x); }
    T operator()(Supports_Condition_Ptr x)     { return static_cast<D*>(this)->fallback(x); }
    T operator()(Supports_Operator_Ptr x)      { return static_cast<D*>(this)->fallback(x); }
    T operator()(Supports_Negation_Ptr x)      { return static_cast<D*>(this)->fallback(x); }
    T operator()(Supports_Declaration_Ptr x)   { return static_cast<D*>(this)->fallback(x); }
    T operator()(Supports_Interpolation_Ptr x) { return static_cast<D*>(this)->fallback(x); }
    T operator()(Media_Query_Ptr x)            { return static_cast<D*>(this)->fallback(x); }
    T operator()(Media_Query_Expression_Ptr x) { return static_cast<D*>(this)->fallback(x); }
    T operator()(At_Root_Query_Ptr x)          { return static_cast<D*>(this)->fallback(x); }
    T operator()(Parent_Selector_Ptr x)        { return static_cast<D*>(this)->fallback(x); }
    T operator()(Parent_Reference_Ptr x)        { return static_cast<D*>(this)->fallback(x); }
    // parameters and arguments
    T operator()(Parameter_Ptr x)              { return static_cast<D*>(this)->fallback(x); }
    T operator()(Parameters_Ptr x)             { return static_cast<D*>(this)->fallback(x); }
    T operator()(Argument_Ptr x)               { return static_cast<D*>(this)->fallback(x); }
    T operator()(Arguments_Ptr x)              { return static_cast<D*>(this)->fallback(x); }
    // selectors
    T operator()(Selector_Schema_Ptr x)        { return static_cast<D*>(this)->fallback(x); }
    T operator()(Placeholder_Selector_Ptr x)   { return static_cast<D*>(this)->fallback(x); }
    T operator()(Type_Selector_Ptr x)       { return static_cast<D*>(this)->fallback(x); }
    T operator()(Class_Selector_Ptr x)         { return static_cast<D*>(this)->fallback(x); }
    T operator()(Id_Selector_Ptr x)            { return static_cast<D*>(this)->fallback(x); }
    T operator()(Attribute_Selector_Ptr x)     { return static_cast<D*>(this)->fallback(x); }
    T operator()(Pseudo_Selector_Ptr x)        { return static_cast<D*>(this)->fallback(x); }
    T operator()(Wrapped_Selector_Ptr x)       { return static_cast<D*>(this)->fallback(x); }
    T operator()(Compound_Selector_Ptr x){ return static_cast<D*>(this)->fallback(x); }
    T operator()(Complex_Selector_Ptr x)      { return static_cast<D*>(this)->fallback(x); }
    T operator()(Selector_List_Ptr x) { return static_cast<D*>(this)->fallback(x); }

    // fallback with specific type U
    // will be called if not overloaded
    template <typename U> T fallback(U x)
    {
      throw std::runtime_error(
        std::string(typeid(*this).name()) + ": CRTP not implemented for " + typeid(x).name());
    }

  };

}

#endif
#ifndef SASS_EMITTER_H
#define SASS_EMITTER_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {
  class Context;

  class Emitter {

    public:
      Emitter(struct Sass_Output_Options& opt);
      virtual ~Emitter() { }

    protected:
      OutputBuffer wbuf;
    public:
      const std::string& buffer(void) { return wbuf.buffer; }
      const SourceMap smap(void) { return wbuf.smap; }
      const OutputBuffer output(void) { return wbuf; }
      // proxy methods for source maps
      void add_source_index(size_t idx);
      void set_filename(const std::string& str);
      void add_open_mapping(AST_Node_Ptr_Const node);
      void add_close_mapping(AST_Node_Ptr_Const node);
      void schedule_mapping(AST_Node_Ptr_Const node);
      std::string render_srcmap(Context &ctx);
      ParserState remap(const ParserState& pstate);

    public:
      struct Sass_Output_Options& opt;
      size_t indentation;
      size_t scheduled_space;
      size_t scheduled_linefeed;
      bool scheduled_delimiter;
      AST_Node_Ptr_Const scheduled_crutch;
      AST_Node_Ptr_Const scheduled_mapping;

    public:
      // output strings different in custom css properties
      bool in_custom_property;
      // output strings different in comments
      bool in_comment;
      // selector list does not get linefeeds
      bool in_wrapped;
      // lists always get a space after delimiter
      bool in_media_block;
      // nested list must not have parentheses
      bool in_declaration;
      // nested lists need parentheses
      bool in_space_array;
      bool in_comma_array;

    public:
      // return buffer as std::string
      std::string get_buffer(void);
      // flush scheduled space/linefeed
      Sass_Output_Style output_style(void) const;
      // add outstanding linefeed
      void finalize(bool final = true);
      // flush scheduled space/linefeed
      void flush_schedules(void);
      // prepend some text or token to the buffer
      void prepend_string(const std::string& text);
      void prepend_output(const OutputBuffer& out);
      // append some text or token to the buffer
      void append_string(const std::string& text);
      // append a single character to buffer
      void append_char(const char chr);
      // append some white-space only text
      void append_wspace(const std::string& text);
      // append some text or token to the buffer
      // this adds source-mappings for node start and end
      void append_token(const std::string& text, const AST_Node_Ptr node);
      // query last appended character
      char last_char();

    public: // syntax sugar
      void append_indentation();
      void append_optional_space(void);
      void append_mandatory_space(void);
      void append_special_linefeed(void);
      void append_optional_linefeed(void);
      void append_mandatory_linefeed(void);
      void append_scope_opener(AST_Node_Ptr node = 0);
      void append_scope_closer(AST_Node_Ptr node = 0);
      void append_comma_separator(void);
      void append_colon_separator(void);
      void append_delimiter(void);

  };

}

#endif

namespace Sass {
  class Context;

  class Inspect : public Operation_CRTP<void, Inspect>, public Emitter {
  protected:
    // import all the class-specific methods and override as desired
    using Operation_CRTP<void, Inspect>::operator();

  public:

    Inspect(const Emitter& emi);
    virtual ~Inspect();

    // statements
    virtual void operator()(Block_Ptr);
    virtual void operator()(Ruleset_Ptr);
    virtual void operator()(Bubble_Ptr);
    virtual void operator()(Supports_Block_Ptr);
    virtual void operator()(Media_Block_Ptr);
    virtual void operator()(At_Root_Block_Ptr);
    virtual void operator()(Directive_Ptr);
    virtual void operator()(Keyframe_Rule_Ptr);
    virtual void operator()(Declaration_Ptr);
    virtual void operator()(Assignment_Ptr);
    virtual void operator()(Import_Ptr);
    virtual void operator()(Import_Stub_Ptr);
    virtual void operator()(Warning_Ptr);
    virtual void operator()(Error_Ptr);
    virtual void operator()(Debug_Ptr);
    virtual void operator()(Comment_Ptr);
    virtual void operator()(If_Ptr);
    virtual void operator()(For_Ptr);
    virtual void operator()(Each_Ptr);
    virtual void operator()(While_Ptr);
    virtual void operator()(Return_Ptr);
    virtual void operator()(Extension_Ptr);
    virtual void operator()(Definition_Ptr);
    virtual void operator()(Mixin_Call_Ptr);
    virtual void operator()(Content_Ptr);
    // expressions
    virtual void operator()(Map_Ptr);
    virtual void operator()(Function_Ptr);
    virtual void operator()(List_Ptr);
    virtual void operator()(Binary_Expression_Ptr);
    virtual void operator()(Unary_Expression_Ptr);
    virtual void operator()(Function_Call_Ptr);
    // virtual void operator()(Custom_Warning_Ptr);
    // virtual void operator()(Custom_Error_Ptr);
    virtual void operator()(Variable_Ptr);
    virtual void operator()(Number_Ptr);
    virtual void operator()(Color_RGBA_Ptr);
    virtual void operator()(Color_HSLA_Ptr);
    virtual void operator()(Boolean_Ptr);
    virtual void operator()(String_Schema_Ptr);
    virtual void operator()(String_Constant_Ptr);
    virtual void operator()(String_Quoted_Ptr);
    virtual void operator()(Custom_Error_Ptr);
    virtual void operator()(Custom_Warning_Ptr);
    virtual void operator()(Supports_Operator_Ptr);
    virtual void operator()(Supports_Negation_Ptr);
    virtual void operator()(Supports_Declaration_Ptr);
    virtual void operator()(Supports_Interpolation_Ptr);
    virtual void operator()(Media_Query_Ptr);
    virtual void operator()(Media_Query_Expression_Ptr);
    virtual void operator()(At_Root_Query_Ptr);
    virtual void operator()(Null_Ptr);
    virtual void operator()(Parent_Selector_Ptr p);
    // parameters and arguments
    virtual void operator()(Parameter_Ptr);
    virtual void operator()(Parameters_Ptr);
    virtual void operator()(Argument_Ptr);
    virtual void operator()(Arguments_Ptr);
    // selectors
    virtual void operator()(Selector_Schema_Ptr);
    virtual void operator()(Placeholder_Selector_Ptr);
    virtual void operator()(Type_Selector_Ptr);
    virtual void operator()(Class_Selector_Ptr);
    virtual void operator()(Id_Selector_Ptr);
    virtual void operator()(Attribute_Selector_Ptr);
    virtual void operator()(Pseudo_Selector_Ptr);
    virtual void operator()(Wrapped_Selector_Ptr);
    virtual void operator()(Compound_Selector_Ptr);
    virtual void operator()(Complex_Selector_Ptr);
    virtual void operator()(Selector_List_Ptr);

    virtual std::string lbracket(List_Ptr);
    virtual std::string rbracket(List_Ptr);

  };

}
#endif

namespace Sass {
  class Context;

  class Output : public Inspect {
  protected:
    using Inspect::operator();

  public:
    Output(Sass_Output_Options& opt);
    virtual ~Output();

  protected:
    std::string charset;
    std::vector<AST_Node_Ptr> top_nodes;

  public:
    OutputBuffer get_buffer(void);

    virtual void operator()(Map_Ptr);
    virtual void operator()(Ruleset_Ptr);
    virtual void operator()(Supports_Block_Ptr);
    virtual void operator()(Media_Block_Ptr);
    virtual void operator()(Directive_Ptr);
    virtual void operator()(Keyframe_Rule_Ptr);
    virtual void operator()(Import_Ptr);
    virtual void operator()(Comment_Ptr);
    virtual void operator()(Number_Ptr);
    virtual void operator()(String_Quoted_Ptr);
    virtual void operator()(String_Constant_Ptr);

    void fallback_impl(AST_Node_Ptr n);

  };

}

#endif
#ifndef SASS_PLUGINS_H
#define SASS_PLUGINS_H

#ifndef SASS_UTF8_STRING_H
#define SASS_UTF8_STRING_H

// Copyright 2006 Nemanja Trifunovic

/*
Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/


#ifndef UTF8_FOR_CPP_2675DCD0_9480_4c0c_B92A_CC14C027B731
#define UTF8_FOR_CPP_2675DCD0_9480_4c0c_B92A_CC14C027B731

// Copyright 2006-2016 Nemanja Trifunovic

/*
Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/


#ifndef UTF8_FOR_CPP_CHECKED_H_2675DCD0_9480_4c0c_B92A_CC14C027B731
#define UTF8_FOR_CPP_CHECKED_H_2675DCD0_9480_4c0c_B92A_CC14C027B731

// Copyright 2006 Nemanja Trifunovic

/*
Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/


#ifndef UTF8_FOR_CPP_CORE_H_2675DCD0_9480_4c0c_B92A_CC14C027B731
#define UTF8_FOR_CPP_CORE_H_2675DCD0_9480_4c0c_B92A_CC14C027B731


namespace utf8
{
    // The typedefs for 8-bit, 16-bit and 32-bit unsigned integers
    // You may need to change them to match your system.
    // These typedefs have the same names as ones from cstdint, or boost/cstdint
    typedef unsigned char   uint8_t;
    typedef unsigned short  uint16_t;
    typedef unsigned int    uint32_t;

// Helper code - not intended to be directly called by the library users. May be changed at any time
namespace internal
{
    // Unicode constants
    // Leading (high) surrogates: 0xd800 - 0xdbff
    // Trailing (low) surrogates: 0xdc00 - 0xdfff
    const uint16_t LEAD_SURROGATE_MIN  = 0xd800u;
    const uint16_t LEAD_SURROGATE_MAX  = 0xdbffu;
    const uint16_t TRAIL_SURROGATE_MIN = 0xdc00u;
    const uint16_t TRAIL_SURROGATE_MAX = 0xdfffu;
    const uint16_t LEAD_OFFSET         = LEAD_SURROGATE_MIN - (0x10000 >> 10);
    const uint32_t SURROGATE_OFFSET    = 0x10000u - (LEAD_SURROGATE_MIN << 10) - TRAIL_SURROGATE_MIN;

    // Maximum valid value for a Unicode code point
    const uint32_t CODE_POINT_MAX      = 0x0010ffffu;

    template<typename octet_type>
    inline uint8_t mask8(octet_type oc)
    {
        return static_cast<uint8_t>(0xff & oc);
    }
    template<typename u16_type>
    inline uint16_t mask16(u16_type oc)
    {
        return static_cast<uint16_t>(0xffff & oc);
    }
    template<typename octet_type>
    inline bool is_trail(octet_type oc)
    {
        return ((utf8::internal::mask8(oc) >> 6) == 0x2);
    }

    template <typename u16>
    inline bool is_lead_surrogate(u16 cp)
    {
        return (cp >= LEAD_SURROGATE_MIN && cp <= LEAD_SURROGATE_MAX);
    }

    template <typename u16>
    inline bool is_trail_surrogate(u16 cp)
    {
        return (cp >= TRAIL_SURROGATE_MIN && cp <= TRAIL_SURROGATE_MAX);
    }

    template <typename u16>
    inline bool is_surrogate(u16 cp)
    {
        return (cp >= LEAD_SURROGATE_MIN && cp <= TRAIL_SURROGATE_MAX);
    }

    template <typename u32>
    inline bool is_code_point_valid(u32 cp)
    {
        return (cp <= CODE_POINT_MAX && !utf8::internal::is_surrogate(cp));
    }

    template <typename octet_iterator>
    inline typename std::iterator_traits<octet_iterator>::difference_type
    sequence_length(octet_iterator lead_it)
    {
        uint8_t lead = utf8::internal::mask8(*lead_it);
        if (lead < 0x80)
            return 1;
        else if ((lead >> 5) == 0x6)
            return 2;
        else if ((lead >> 4) == 0xe)
            return 3;
        else if ((lead >> 3) == 0x1e)
            return 4;
        else
            return 0;
    }

    template <typename octet_difference_type>
    inline bool is_overlong_sequence(uint32_t cp, octet_difference_type length)
    {
        if (cp < 0x80) {
            if (length != 1)
                return true;
        }
        else if (cp < 0x800) {
            if (length != 2)
                return true;
        }
        else if (cp < 0x10000) {
            if (length != 3)
                return true;
        }

        return false;
    }

    enum utf_error {UTF8_OK, NOT_ENOUGH_ROOM, INVALID_LEAD, INCOMPLETE_SEQUENCE, OVERLONG_SEQUENCE, INVALID_CODE_POINT};

    /// Helper for get_sequence_x
    template <typename octet_iterator>
    utf_error increase_safely(octet_iterator& it, octet_iterator end)
    {
        if (++it == end)
            return NOT_ENOUGH_ROOM;

        if (!utf8::internal::is_trail(*it))
            return INCOMPLETE_SEQUENCE;

        return UTF8_OK;
    }

    #define UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(IT, END) {utf_error ret = increase_safely(IT, END); if (ret != UTF8_OK) return ret;}

    /// get_sequence_x functions decode utf-8 sequences of the length x
    template <typename octet_iterator>
    utf_error get_sequence_1(octet_iterator& it, octet_iterator end, uint32_t& code_point)
    {
        if (it == end)
            return NOT_ENOUGH_ROOM;

        code_point = utf8::internal::mask8(*it);

        return UTF8_OK;
    }

    template <typename octet_iterator>
    utf_error get_sequence_2(octet_iterator& it, octet_iterator end, uint32_t& code_point)
    {
        if (it == end)
            return NOT_ENOUGH_ROOM;

        code_point = utf8::internal::mask8(*it);

        UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(it, end)

        code_point = ((code_point << 6) & 0x7ff) + ((*it) & 0x3f);

        return UTF8_OK;
    }

    template <typename octet_iterator>
    utf_error get_sequence_3(octet_iterator& it, octet_iterator end, uint32_t& code_point)
    {
        if (it == end)
            return NOT_ENOUGH_ROOM;

        code_point = utf8::internal::mask8(*it);

        UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(it, end)

        code_point = ((code_point << 12) & 0xffff) + ((utf8::internal::mask8(*it) << 6) & 0xfff);

        UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(it, end)

        code_point += (*it) & 0x3f;

        return UTF8_OK;
    }

    template <typename octet_iterator>
    utf_error get_sequence_4(octet_iterator& it, octet_iterator end, uint32_t& code_point)
    {
        if (it == end)
           return NOT_ENOUGH_ROOM;

        code_point = utf8::internal::mask8(*it);

        UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(it, end)

        code_point = ((code_point << 18) & 0x1fffff) + ((utf8::internal::mask8(*it) << 12) & 0x3ffff);

        UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(it, end)

        code_point += (utf8::internal::mask8(*it) << 6) & 0xfff;

        UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR(it, end)

        code_point += (*it) & 0x3f;

        return UTF8_OK;
    }

    #undef UTF8_CPP_INCREASE_AND_RETURN_ON_ERROR

    template <typename octet_iterator>
    utf_error validate_next(octet_iterator& it, octet_iterator end, uint32_t& code_point)
    {
        if (it == end)
            return NOT_ENOUGH_ROOM;

        // Save the original value of it so we can go back in case of failure
        // Of course, it does not make much sense with i.e. stream iterators
        octet_iterator original_it = it;

        uint32_t cp = 0;
        // Determine the sequence length based on the lead octet
        typedef typename std::iterator_traits<octet_iterator>::difference_type octet_difference_type;
        const octet_difference_type length = utf8::internal::sequence_length(it);

        // Get trail octets and calculate the code point
        utf_error err = UTF8_OK;
        switch (length) {
            case 0:
                return INVALID_LEAD;
            case 1:
                err = utf8::internal::get_sequence_1(it, end, cp);
                break;
            case 2:
                err = utf8::internal::get_sequence_2(it, end, cp);
            break;
            case 3:
                err = utf8::internal::get_sequence_3(it, end, cp);
            break;
            case 4:
                err = utf8::internal::get_sequence_4(it, end, cp);
            break;
        }

        if (err == UTF8_OK) {
            // Decoding succeeded. Now, security checks...
            if (utf8::internal::is_code_point_valid(cp)) {
                if (!utf8::internal::is_overlong_sequence(cp, length)){
                    // Passed! Return here.
                    code_point = cp;
                    ++it;
                    return UTF8_OK;
                }
                else
                    err = OVERLONG_SEQUENCE;
            }
            else
                err = INVALID_CODE_POINT;
        }

        // Failure branch - restore the original value of the iterator
        it = original_it;
        return err;
    }

    template <typename octet_iterator>
    inline utf_error validate_next(octet_iterator& it, octet_iterator end) {
        uint32_t ignored;
        return utf8::internal::validate_next(it, end, ignored);
    }

} // namespace internal

    /// The library API - functions intended to be called by the users

    // Byte order mark
    const uint8_t bom[] = {0xef, 0xbb, 0xbf};

    template <typename octet_iterator>
    octet_iterator find_invalid(octet_iterator start, octet_iterator end)
    {
        octet_iterator result = start;
        while (result != end) {
            utf8::internal::utf_error err_code = utf8::internal::validate_next(result, end);
            if (err_code != internal::UTF8_OK)
                return result;
        }
        return result;
    }

    template <typename octet_iterator>
    inline bool is_valid(octet_iterator start, octet_iterator end)
    {
        return (utf8::find_invalid(start, end) == end);
    }

    template <typename octet_iterator>
    inline bool starts_with_bom (octet_iterator it, octet_iterator end)
    {
        return (
            ((it != end) && (utf8::internal::mask8(*it++)) == bom[0]) &&
            ((it != end) && (utf8::internal::mask8(*it++)) == bom[1]) &&
            ((it != end) && (utf8::internal::mask8(*it))   == bom[2])
           );
    }

    //Deprecated in release 2.3
    template <typename octet_iterator>
    inline bool is_bom (octet_iterator it)
    {
        return (
            (utf8::internal::mask8(*it++)) == bom[0] &&
            (utf8::internal::mask8(*it++)) == bom[1] &&
            (utf8::internal::mask8(*it))   == bom[2]
           );
    }
} // namespace utf8

#endif // header guard



namespace utf8
{
    // Base for the exceptions that may be thrown from the library
    class exception : public ::std::exception {
    };

    // Exceptions that may be thrown from the library functions.
    class invalid_code_point : public exception {
        uint32_t cp;
    public:
        invalid_code_point(uint32_t codepoint) : cp(codepoint) {}
        virtual const char* what() const throw() { return "Invalid code point"; }
        uint32_t code_point() const {return cp;}
    };

    class invalid_utf8 : public exception {
        uint8_t u8;
    public:
        invalid_utf8 (uint8_t u) : u8(u) {}
        virtual const char* what() const throw() { return "Invalid UTF-8"; }
        uint8_t utf8_octet() const {return u8;}
    };

    class invalid_utf16 : public exception {
        uint16_t u16;
    public:
        invalid_utf16 (uint16_t u) : u16(u) {}
        virtual const char* what() const throw() { return "Invalid UTF-16"; }
        uint16_t utf16_word() const {return u16;}
    };

    class not_enough_room : public exception {
    public:
        virtual const char* what() const throw() { return "Not enough space"; }
    };

    /// The library API - functions intended to be called by the users

    template <typename octet_iterator>
    octet_iterator append(uint32_t cp, octet_iterator result)
    {
        if (!utf8::internal::is_code_point_valid(cp))
            throw invalid_code_point(cp);

        if (cp < 0x80)                        // one octet
            *(result++) = static_cast<uint8_t>(cp);
        else if (cp < 0x800) {                // two octets
            *(result++) = static_cast<uint8_t>((cp >> 6)            | 0xc0);
            *(result++) = static_cast<uint8_t>((cp & 0x3f)          | 0x80);
        }
        else if (cp < 0x10000) {              // three octets
            *(result++) = static_cast<uint8_t>((cp >> 12)           | 0xe0);
            *(result++) = static_cast<uint8_t>(((cp >> 6) & 0x3f)   | 0x80);
            *(result++) = static_cast<uint8_t>((cp & 0x3f)          | 0x80);
        }
        else {                                // four octets
            *(result++) = static_cast<uint8_t>((cp >> 18)           | 0xf0);
            *(result++) = static_cast<uint8_t>(((cp >> 12) & 0x3f)  | 0x80);
            *(result++) = static_cast<uint8_t>(((cp >> 6) & 0x3f)   | 0x80);
            *(result++) = static_cast<uint8_t>((cp & 0x3f)          | 0x80);
        }
        return result;
    }

    template <typename octet_iterator, typename output_iterator>
    output_iterator replace_invalid(octet_iterator start, octet_iterator end, output_iterator out, uint32_t replacement)
    {
        while (start != end) {
            octet_iterator sequence_start = start;
            internal::utf_error err_code = utf8::internal::validate_next(start, end);
            switch (err_code) {
                case internal::UTF8_OK :
                    for (octet_iterator it = sequence_start; it != start; ++it)
                        *out++ = *it;
                    break;
                case internal::NOT_ENOUGH_ROOM:
                    out = utf8::append (replacement, out);
                    start = end;
                    break;
                case internal::INVALID_LEAD:
                    out = utf8::append (replacement, out);
                    ++start;
                    break;
                case internal::INCOMPLETE_SEQUENCE:
                case internal::OVERLONG_SEQUENCE:
                case internal::INVALID_CODE_POINT:
                    out = utf8::append (replacement, out);
                    ++start;
                    // just one replacement mark for the sequence
                    while (start != end && utf8::internal::is_trail(*start))
                        ++start;
                    break;
            }
        }
        return out;
    }

    template <typename octet_iterator, typename output_iterator>
    inline output_iterator replace_invalid(octet_iterator start, octet_iterator end, output_iterator out)
    {
        static const uint32_t replacement_marker = utf8::internal::mask16(0xfffd);
        return utf8::replace_invalid(start, end, out, replacement_marker);
    }

    template <typename octet_iterator>
    uint32_t next(octet_iterator& it, octet_iterator end)
    {
        uint32_t cp = 0;
        internal::utf_error err_code = utf8::internal::validate_next(it, end, cp);
        switch (err_code) {
            case internal::UTF8_OK :
                break;
            case internal::NOT_ENOUGH_ROOM :
                throw not_enough_room();
            case internal::INVALID_LEAD :
            case internal::INCOMPLETE_SEQUENCE :
            case internal::OVERLONG_SEQUENCE :
                throw invalid_utf8(*it);
            case internal::INVALID_CODE_POINT :
                throw invalid_code_point(cp);
        }
        return cp;
    }

    template <typename octet_iterator>
    uint32_t peek_next(octet_iterator it, octet_iterator end)
    {
        return utf8::next(it, end);
    }

    template <typename octet_iterator>
    uint32_t prior(octet_iterator& it, octet_iterator start)
    {
        // can't do much if it == start
        if (it == start)
            throw not_enough_room();

        octet_iterator end = it;
        // Go back until we hit either a lead octet or start
        while (utf8::internal::is_trail(*(--it)))
            if (it == start)
                throw invalid_utf8(*it); // error - no lead byte in the sequence
        return utf8::peek_next(it, end);
    }

    /// Deprecated in versions that include "prior"
    template <typename octet_iterator>
    uint32_t previous(octet_iterator& it, octet_iterator pass_start)
    {
        octet_iterator end = it;
        while (utf8::internal::is_trail(*(--it)))
            if (it == pass_start)
                throw invalid_utf8(*it); // error - no lead byte in the sequence
        octet_iterator temp = it;
        return utf8::next(temp, end);
    }

    template <typename octet_iterator, typename distance_type>
    void advance (octet_iterator& it, distance_type n, octet_iterator end)
    {
        for (distance_type i = 0; i < n; ++i)
            utf8::next(it, end);
    }

    template <typename octet_iterator, typename distance_type>
    void retreat (octet_iterator& it, distance_type n, octet_iterator end)
    {
        for (distance_type i = 0; i < n; ++i)
            utf8::prior(it, end);
    }

    template <typename octet_iterator>
    typename std::iterator_traits<octet_iterator>::difference_type
    distance (octet_iterator first, octet_iterator last)
    {
        typename std::iterator_traits<octet_iterator>::difference_type dist;
        for (dist = 0; first < last; ++dist)
            utf8::next(first, last);
        return dist;
    }

    template <typename u16bit_iterator, typename octet_iterator>
    octet_iterator utf16to8 (u16bit_iterator start, u16bit_iterator end, octet_iterator result)
    {
        while (start != end) {
            uint32_t cp = utf8::internal::mask16(*start++);
            // Take care of surrogate pairs first
            if (utf8::internal::is_lead_surrogate(cp)) {
                if (start != end) {
                    uint32_t trail_surrogate = utf8::internal::mask16(*start++);
                    if (utf8::internal::is_trail_surrogate(trail_surrogate))
                        cp = (cp << 10) + trail_surrogate + internal::SURROGATE_OFFSET;
                    else
                        throw invalid_utf16(static_cast<uint16_t>(trail_surrogate));
                }
                else
                    throw invalid_utf16(static_cast<uint16_t>(cp));

            }
            // Lone trail surrogate
            else if (utf8::internal::is_trail_surrogate(cp))
                throw invalid_utf16(static_cast<uint16_t>(cp));

            result = utf8::append(cp, result);
        }
        return result;
    }

    template <typename u16bit_iterator, typename octet_iterator>
    u16bit_iterator utf8to16 (octet_iterator start, octet_iterator end, u16bit_iterator result)
    {
        while (start < end) {
            uint32_t cp = utf8::next(start, end);
            if (cp > 0xffff) { //make a surrogate pair
                *result++ = static_cast<uint16_t>((cp >> 10)   + internal::LEAD_OFFSET);
                *result++ = static_cast<uint16_t>((cp & 0x3ff) + internal::TRAIL_SURROGATE_MIN);
            }
            else
                *result++ = static_cast<uint16_t>(cp);
        }
        return result;
    }

    template <typename octet_iterator, typename u32bit_iterator>
    octet_iterator utf32to8 (u32bit_iterator start, u32bit_iterator end, octet_iterator result)
    {
        while (start != end)
            result = utf8::append(*(start++), result);

        return result;
    }

    template <typename octet_iterator, typename u32bit_iterator>
    u32bit_iterator utf8to32 (octet_iterator start, octet_iterator end, u32bit_iterator result)
    {
        while (start < end)
            (*result++) = utf8::next(start, end);

        return result;
    }

    // The iterator class
    template <typename octet_iterator>
    class iterator : public std::iterator <std::bidirectional_iterator_tag, uint32_t> {
      octet_iterator it;
      octet_iterator range_start;
      octet_iterator range_end;
      public:
      iterator () {}
      explicit iterator (const octet_iterator& octet_it,
                         const octet_iterator& rangestart,
                         const octet_iterator& rangeend) :
               it(octet_it), range_start(rangestart), range_end(rangeend)
      {
          if (it < range_start || it > range_end)
              throw std::out_of_range("Invalid utf-8 iterator position");
      }
      // the default "big three" are OK
      octet_iterator base () const { return it; }
      uint32_t operator * () const
      {
          octet_iterator temp = it;
          return utf8::next(temp, range_end);
      }
      bool operator == (const iterator& rhs) const
      {
          if (range_start != rhs.range_start || range_end != rhs.range_end)
              throw std::logic_error("Comparing utf-8 iterators defined with different ranges");
          return (it == rhs.it);
      }
      bool operator != (const iterator& rhs) const
      {
          return !(operator == (rhs));
      }
      iterator& operator ++ ()
      {
          utf8::next(it, range_end);
          return *this;
      }
      iterator operator ++ (int)
      {
          iterator temp = *this;
          utf8::next(it, range_end);
          return temp;
      }
      iterator& operator -- ()
      {
          utf8::prior(it, range_start);
          return *this;
      }
      iterator operator -- (int)
      {
          iterator temp = *this;
          utf8::prior(it, range_start);
          return temp;
      }
    }; // class iterator

} // namespace utf8

#endif //header guard


// Copyright 2006 Nemanja Trifunovic

/*
Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
*/


#ifndef UTF8_FOR_CPP_UNCHECKED_H_2675DCD0_9480_4c0c_B92A_CC14C027B731
#define UTF8_FOR_CPP_UNCHECKED_H_2675DCD0_9480_4c0c_B92A_CC14C027B731


namespace utf8
{
    namespace unchecked
    {
        template <typename octet_iterator>
        octet_iterator append(uint32_t cp, octet_iterator result)
        {
            if (cp < 0x80)                        // one octet
                *(result++) = static_cast<uint8_t>(cp);
            else if (cp < 0x800) {                // two octets
                *(result++) = static_cast<uint8_t>((cp >> 6)          | 0xc0);
                *(result++) = static_cast<uint8_t>((cp & 0x3f)        | 0x80);
            }
            else if (cp < 0x10000) {              // three octets
                *(result++) = static_cast<uint8_t>((cp >> 12)         | 0xe0);
                *(result++) = static_cast<uint8_t>(((cp >> 6) & 0x3f) | 0x80);
                *(result++) = static_cast<uint8_t>((cp & 0x3f)        | 0x80);
            }
            else {                                // four octets
                *(result++) = static_cast<uint8_t>((cp >> 18)         | 0xf0);
                *(result++) = static_cast<uint8_t>(((cp >> 12) & 0x3f)| 0x80);
                *(result++) = static_cast<uint8_t>(((cp >> 6) & 0x3f) | 0x80);
                *(result++) = static_cast<uint8_t>((cp & 0x3f)        | 0x80);
            }
            return result;
        }

        template <typename octet_iterator>
        uint32_t next(octet_iterator& it)
        {
            uint32_t cp = utf8::internal::mask8(*it);
            typename std::iterator_traits<octet_iterator>::difference_type length = utf8::internal::sequence_length(it);
            switch (length) {
                case 1:
                    break;
                case 2:
                    it++;
                    cp = ((cp << 6) & 0x7ff) + ((*it) & 0x3f);
                    break;
                case 3:
                    ++it;
                    cp = ((cp << 12) & 0xffff) + ((utf8::internal::mask8(*it) << 6) & 0xfff);
                    ++it;
                    cp += (*it) & 0x3f;
                    break;
                case 4:
                    ++it;
                    cp = ((cp << 18) & 0x1fffff) + ((utf8::internal::mask8(*it) << 12) & 0x3ffff);
                    ++it;
                    cp += (utf8::internal::mask8(*it) << 6) & 0xfff;
                    ++it;
                    cp += (*it) & 0x3f;
                    break;
            }
            ++it;
            return cp;
        }

        template <typename octet_iterator>
        uint32_t peek_next(octet_iterator it)
        {
            return utf8::unchecked::next(it);
        }

        template <typename octet_iterator>
        uint32_t prior(octet_iterator& it)
        {
            while (utf8::internal::is_trail(*(--it))) ;
            octet_iterator temp = it;
            return utf8::unchecked::next(temp);
        }

        // Deprecated in versions that include prior, but only for the sake of consistency (see utf8::previous)
        template <typename octet_iterator>
        inline uint32_t previous(octet_iterator& it)
        {
            return utf8::unchecked::prior(it);
        }

        template <typename octet_iterator, typename distance_type>
        void advance (octet_iterator& it, distance_type n)
        {
            for (distance_type i = 0; i < n; ++i)
                utf8::unchecked::next(it);
        }

        template <typename octet_iterator, typename distance_type>
        void retreat (octet_iterator& it, distance_type n)
        {
            for (distance_type i = 0; i < n; ++i)
                utf8::unchecked::prior(it);
        }

        template <typename octet_iterator>
        typename std::iterator_traits<octet_iterator>::difference_type
        distance (octet_iterator first, octet_iterator last)
        {
            typename std::iterator_traits<octet_iterator>::difference_type dist;
            for (dist = 0; first < last; ++dist)
                utf8::unchecked::next(first);
            return dist;
        }

        template <typename u16bit_iterator, typename octet_iterator>
        octet_iterator utf16to8 (u16bit_iterator start, u16bit_iterator end, octet_iterator result)
        {
            while (start != end) {
                uint32_t cp = utf8::internal::mask16(*start++);
            // Take care of surrogate pairs first
                if (utf8::internal::is_lead_surrogate(cp)) {
                    uint32_t trail_surrogate = utf8::internal::mask16(*start++);
                    cp = (cp << 10) + trail_surrogate + internal::SURROGATE_OFFSET;
                }
                result = utf8::unchecked::append(cp, result);
            }
            return result;
        }

        template <typename u16bit_iterator, typename octet_iterator>
        u16bit_iterator utf8to16 (octet_iterator start, octet_iterator end, u16bit_iterator result)
        {
            while (start < end) {
                uint32_t cp = utf8::unchecked::next(start);
                if (cp > 0xffff) { //make a surrogate pair
                    *result++ = static_cast<uint16_t>((cp >> 10)   + internal::LEAD_OFFSET);
                    *result++ = static_cast<uint16_t>((cp & 0x3ff) + internal::TRAIL_SURROGATE_MIN);
                }
                else
                    *result++ = static_cast<uint16_t>(cp);
            }
            return result;
        }

        template <typename octet_iterator, typename u32bit_iterator>
        octet_iterator utf32to8 (u32bit_iterator start, u32bit_iterator end, octet_iterator result)
        {
            while (start != end)
                result = utf8::unchecked::append(*(start++), result);

            return result;
        }

        template <typename octet_iterator, typename u32bit_iterator>
        u32bit_iterator utf8to32 (octet_iterator start, octet_iterator end, u32bit_iterator result)
        {
            while (start < end)
                (*result++) = utf8::unchecked::next(start);

            return result;
        }

        // The iterator class
        template <typename octet_iterator>
          class iterator : public std::iterator <std::bidirectional_iterator_tag, uint32_t> {
            octet_iterator it;
            public:
            iterator () {}
            explicit iterator (const octet_iterator& octet_it): it(octet_it) {}
            // the default "big three" are OK
            octet_iterator base () const { return it; }
            uint32_t operator * () const
            {
                octet_iterator temp = it;
                return utf8::unchecked::next(temp);
            }
            bool operator == (const iterator& rhs) const
            {
                return (it == rhs.it);
            }
            bool operator != (const iterator& rhs) const
            {
                return !(operator == (rhs));
            }
            iterator& operator ++ ()
            {
                ::std::advance(it, utf8::internal::sequence_length(it));
                return *this;
            }
            iterator operator ++ (int)
            {
                iterator temp = *this;
                ::std::advance(it, utf8::internal::sequence_length(it));
                return temp;
            }
            iterator& operator -- ()
            {
                utf8::unchecked::prior(it);
                return *this;
            }
            iterator operator -- (int)
            {
                iterator temp = *this;
                utf8::unchecked::prior(it);
                return temp;
            }
          }; // class iterator

    } // namespace utf8::unchecked
} // namespace utf8


#endif // header guard


#endif // header guard

namespace Sass {
  namespace UTF_8 {

    // naming conventions:
    // offset: raw byte offset (0 based)
    // position: code point offset (0 based)
    // index: code point offset (1 based or negative)

    // function that will count the number of code points (utf-8 characters) from the beginning to the given end
    size_t code_point_count(const std::string& str, size_t start, size_t end);
    size_t code_point_count(const std::string& str);

    // function that will return the byte offset of a code point in a
    size_t offset_at_position(const std::string& str, size_t position);

    // function that returns number of bytes in a character in a string
    size_t code_point_size_at_offset(const std::string& str, size_t offset);

    // function that will return a normalized index, given a crazy one
    size_t normalize_index(int index, size_t len);

    #ifdef _WIN32
    // functions to handle unicode paths on windows
    std::string convert_from_utf16(const std::wstring& wstr);
    std::wstring convert_to_utf16(const std::string& str);
    #endif

  }
}

#endif

#ifdef _WIN32

  #define LOAD_LIB(var, path) HMODULE var = LoadLibraryW(UTF_8::convert_to_utf16(path).c_str())
  #define LOAD_LIB_WCHR(var, path_wide_str) HMODULE var = LoadLibraryW(path_wide_str.c_str())
  #define LOAD_LIB_FN(type, var, name) type var = (type) GetProcAddress(plugin, name)
  #define CLOSE_LIB(var) FreeLibrary(var)

  #ifndef dlerror
  #define dlerror() 0
  #endif

#else

  #define LOAD_LIB(var, path) void* var = dlopen(path.c_str(), RTLD_LAZY)
  #define LOAD_LIB_FN(type, var, name) type var = (type) dlsym(plugin, name)
  #define CLOSE_LIB(var) dlclose(var)

#endif

namespace Sass {


  class Plugins {

    public: // c-tor
      Plugins(void);
      ~Plugins(void);

    public: // methods
      // load one specific plugin
      bool load_plugin(const std::string& path);
      // load all plugins from a directory
      size_t load_plugins(const std::string& path);

    public: // public accessors
      const std::vector<Sass_Importer_Entry> get_headers(void) { return headers; }
      const std::vector<Sass_Importer_Entry> get_importers(void) { return importers; }
      const std::vector<Sass_Function_Entry> get_functions(void) { return functions; }

    private: // private vars
      std::vector<Sass_Importer_Entry> headers;
      std::vector<Sass_Importer_Entry> importers;
      std::vector<Sass_Function_Entry> functions;

  };

}

#endif


struct Sass_Function;

namespace Sass {

  class Context {
  public:
    void import_url (Import_Ptr imp, std::string load_path, const std::string& ctx_path);
    bool call_headers(const std::string& load_path, const char* ctx_path, ParserState& pstate, Import_Ptr imp)
    { return call_loader(load_path, ctx_path, pstate, imp, c_headers, false); };
    bool call_importers(const std::string& load_path, const char* ctx_path, ParserState& pstate, Import_Ptr imp)
    { return call_loader(load_path, ctx_path, pstate, imp, c_importers, true); };

  private:
    bool call_loader(const std::string& load_path, const char* ctx_path, ParserState& pstate, Import_Ptr imp, std::vector<Sass_Importer_Entry> importers, bool only_one = true);

  public:
    const std::string CWD;
    struct Sass_Options& c_options;
    std::string entry_path;
    size_t head_imports;
    Plugins plugins;
    Output emitter;

    // generic ast node garbage container
    // used to avoid possible circular refs
    CallStack ast_gc;
    // resources add under our control
    // these are guaranteed to be freed
    std::vector<char*> strings;
    std::vector<Resource> resources;
    std::map<const std::string, StyleSheet> sheets;
    Subset_Map subset_map;
    ImporterStack import_stack;
    std::vector<Sass_Callee> callee_stack;
    std::vector<Backtrace> traces;

    struct Sass_Compiler* c_compiler;

    // absolute paths to includes
    std::vector<std::string> included_files;
    // relative includes for sourcemap
    std::vector<std::string> srcmap_links;
    // vectors above have same size

    std::vector<std::string> plugin_paths; // relative paths to load plugins
    std::vector<std::string> include_paths; // lookup paths for includes





    void apply_custom_headers(Block_Obj root, const char* path, ParserState pstate);

    std::vector<Sass_Importer_Entry> c_headers;
    std::vector<Sass_Importer_Entry> c_importers;
    std::vector<Sass_Function_Entry> c_functions;

    void add_c_header(Sass_Importer_Entry header);
    void add_c_importer(Sass_Importer_Entry importer);
    void add_c_function(Sass_Function_Entry function);

    const std::string indent; // String to be used for indentation
    const std::string linefeed; // String to be used for line feeds
    const std::string input_path; // for relative paths in src-map
    const std::string output_path; // for relative paths to the output
    const std::string source_map_file; // path to source map file (enables feature)
    const std::string source_map_root; // path for sourceRoot property (pass-through)

    virtual ~Context();
    Context(struct Sass_Context&);
    virtual Block_Obj parse() = 0;
    virtual Block_Obj compile();
    virtual char* render(Block_Obj root);
    virtual char* render_srcmap();

    void register_resource(const Include&, const Resource&);
    void register_resource(const Include&, const Resource&, ParserState&);
    std::vector<Include> find_includes(const Importer& import);
    Include load_import(const Importer&, ParserState pstate);

    Sass_Output_Style output_style() { return c_options.output_style; };
    std::vector<std::string> get_included_files(bool skip = false, size_t headers = 0);

  private:
    void collect_plugin_paths(const char* paths_str);
    void collect_plugin_paths(string_list* paths_array);
    void collect_include_paths(const char* paths_str);
    void collect_include_paths(string_list* paths_array);
    std::string format_embedded_source_map();
    std::string format_source_mapping_url(const std::string& out_path);


    // void register_built_in_functions(Env* env);
    // void register_function(Signature sig, Native_Function f, Env* env);
    // void register_function(Signature sig, Native_Function f, size_t arity, Env* env);
    // void register_overload_stub(std::string name, Env* env);

  public:
    const std::string& cwd() { return CWD; };
  };

  class File_Context : public Context {
  public:
    File_Context(struct Sass_File_Context& ctx)
    : Context(ctx)
    { }
    virtual ~File_Context();
    virtual Block_Obj parse();
  };

  class Data_Context : public Context {
  public:
    char* source_c_str;
    char* srcmap_c_str;
    Data_Context(struct Sass_Data_Context& ctx)
    : Context(ctx)
    {
      source_c_str       = ctx.source_string;
      srcmap_c_str       = ctx.srcmap_string;
      ctx.source_string = 0; // passed away
      ctx.srcmap_string = 0; // passed away
    }
    virtual ~Data_Context();
    virtual Block_Obj parse();
  };

}

#endif
#ifndef SASS_CONSTANTS_H
#define SASS_CONSTANTS_H

namespace Sass {
  namespace Constants {

    // The maximum call stack that can be created
    extern const unsigned long MaxCallStack;

    // https://developer.mozilla.org/en-US/docs/Web/CSS/Specificity
    // The following list of selectors is by increasing specificity:
    extern const unsigned long Specificity_Star;
    extern const unsigned long Specificity_Universal;
    extern const unsigned long Specificity_Element;
    extern const unsigned long Specificity_Base;
    extern const unsigned long Specificity_Class;
    extern const unsigned long Specificity_Attr;
    extern const unsigned long Specificity_Pseudo;
    extern const unsigned long Specificity_ID;

    // Selector unification order;
    extern const int UnificationOrder_Element;
    extern const int UnificationOrder_Id;
    extern const int UnificationOrder_Class;
    extern const int UnificationOrder_Attribute;
    extern const int UnificationOrder_PseudoClass;
    extern const int UnificationOrder_Wrapped;
    extern const int UnificationOrder_PseudoElement;
    extern const int UnificationOrder_Placeholder;

    // sass keywords
    extern const char at_root_kwd[];
    extern const char import_kwd[];
    extern const char mixin_kwd[];
    extern const char function_kwd[];
    extern const char return_kwd[];
    extern const char include_kwd[];
    extern const char content_kwd[];
    extern const char extend_kwd[];
    extern const char if_kwd[];
    extern const char else_kwd[];
    extern const char if_after_else_kwd[];
    extern const char for_kwd[];
    extern const char from_kwd[];
    extern const char to_kwd[];
    extern const char through_kwd[];
    extern const char each_kwd[];
    extern const char in_kwd[];
    extern const char while_kwd[];
    extern const char warn_kwd[];
    extern const char error_kwd[];
    extern const char debug_kwd[];
    extern const char default_kwd[];
    extern const char global_kwd[];
    extern const char null_kwd[];
    extern const char optional_kwd[];
    extern const char with_kwd[];
    extern const char without_kwd[];
    extern const char all_kwd[];
    extern const char rule_kwd[];

    // css standard units
    extern const char em_kwd[];
    extern const char ex_kwd[];
    extern const char px_kwd[];
    extern const char cm_kwd[];
    extern const char mm_kwd[];
    extern const char pt_kwd[];
    extern const char pc_kwd[];
    extern const char deg_kwd[];
    extern const char rad_kwd[];
    extern const char grad_kwd[];
    extern const char turn_kwd[];
    extern const char ms_kwd[];
    extern const char s_kwd[];
    extern const char Hz_kwd[];
    extern const char kHz_kwd[];

    // vendor prefixes
    extern const char vendor_opera_kwd[];
    extern const char vendor_webkit_kwd[];
    extern const char vendor_mozilla_kwd[];
    extern const char vendor_ms_kwd[];
    extern const char vendor_khtml_kwd[];

    // css functions and keywords
    extern const char charset_kwd[];
    extern const char media_kwd[];
    extern const char supports_kwd[];
    extern const char keyframes_kwd[];
    extern const char only_kwd[];
    extern const char rgb_fn_kwd[];
    extern const char url_fn_kwd[];
    extern const char url_kwd[];
    // extern const char url_prefix_fn_kwd[];
    extern const char important_kwd[];
    extern const char pseudo_not_fn_kwd[];
    extern const char even_kwd[];
    extern const char odd_kwd[];
    extern const char progid_kwd[];
    extern const char expression_kwd[];
    extern const char calc_fn_kwd[];

    // char classes for "regular expressions"
    extern const char almost_any_value_class[];

    // css selector keywords
    extern const char sel_deep_kwd[];

    // css attribute-matching operators
    extern const char tilde_equal[];
    extern const char pipe_equal[];
    extern const char caret_equal[];
    extern const char dollar_equal[];
    extern const char star_equal[];

    // relational & logical operators and constants
    extern const char and_kwd[];
    extern const char or_kwd[];
    extern const char not_kwd[];
    extern const char gt[];
    extern const char gte[];
    extern const char lt[];
    extern const char lte[];
    extern const char eq[];
    extern const char neq[];
    extern const char true_kwd[];
    extern const char false_kwd[];

    // definition keywords
    extern const char using_kwd[];

    // miscellaneous punctuation and delimiters
    extern const char percent_str[];
    extern const char empty_str[];
    extern const char slash_slash[];
    extern const char slash_star[];
    extern const char star_slash[];
    extern const char hash_lbrace[];
    extern const char rbrace[];
    extern const char rparen[];
    extern const char sign_chars[];
    extern const char op_chars[];
    extern const char hyphen[];
    extern const char ellipsis[];
    // extern const char url_space_chars[];

    // type names
    extern const char numeric_name[];
    extern const char number_name[];
    extern const char percentage_name[];
    extern const char dimension_name[];
    extern const char string_name[];
    extern const char bool_name[];
    extern const char color_name[];
    extern const char list_name[];
    extern const char map_name[];
    extern const char arglist_name[];

    // constants for uri parsing (RFC 3986 Appendix A.)
    extern const char uri_chars[];
    extern const char real_uri_chars[];

    // some specific constant character classes
    // they must be static to be useable by lexer
    extern const char static_ops[];
    extern const char selector_list_delims[];
    extern const char complex_selector_delims[];
    extern const char selector_combinator_ops[];
    extern const char attribute_compare_modifiers[];
    extern const char selector_lookahead_ops[];

    // byte order marks
    // (taken from http://en.wikipedia.org/wiki/Byte_order_mark)
    extern const unsigned char utf_8_bom[];
    extern const unsigned char utf_16_bom_be[];
    extern const unsigned char utf_16_bom_le[];
    extern const unsigned char utf_32_bom_be[];
    extern const unsigned char utf_32_bom_le[];
    extern const unsigned char utf_7_bom_1[];
    extern const unsigned char utf_7_bom_2[];
    extern const unsigned char utf_7_bom_3[];
    extern const unsigned char utf_7_bom_4[];
    extern const unsigned char utf_7_bom_5[];
    extern const unsigned char utf_1_bom[];
    extern const unsigned char utf_ebcdic_bom[];
    extern const unsigned char scsu_bom[];
    extern const unsigned char bocu_1_bom[];
    extern const unsigned char gb_18030_bom[];

  }
}

#endif
#ifndef SASS_ERROR_HANDLING_H
#define SASS_ERROR_HANDLING_H


namespace Sass {

  struct Backtrace;

  namespace Exception {

    const std::string def_msg = "Invalid sass detected";
    const std::string def_op_msg = "Undefined operation";
    const std::string def_op_null_msg = "Invalid null operation";
    const std::string def_nesting_limit = "Code too deeply neested";

    class Base : public std::runtime_error {
      protected:
        std::string msg;
        std::string prefix;
      public:
        ParserState pstate;
        Backtraces traces;
      public:
        Base(ParserState pstate, std::string msg, Backtraces traces);
        virtual const char* errtype() const { return prefix.c_str(); }
        virtual const char* what() const throw() { return msg.c_str(); }
        virtual ~Base() throw() {};
    };

    class InvalidSass : public Base {
      public:
        InvalidSass(InvalidSass& other) : Base(other), owned_src(other.owned_src) {
          // Assumes that `this` will outlive `other`.
          other.owned_src = nullptr;
        }

        // Required because the copy constructor's argument is not const.
        // Can't use `std::move` here because we build on Visual Studio 2013.
        InvalidSass(InvalidSass &&other) : Base(other), owned_src(other.owned_src) {
          other.owned_src = nullptr;
        }

        InvalidSass(ParserState pstate, Backtraces traces, std::string msg, char* owned_src = nullptr);
        virtual ~InvalidSass() throw() { sass_free_memory(owned_src); };
        char *owned_src;
    };

    class InvalidParent : public Base {
      protected:
        Selector_Ptr parent;
        Selector_Ptr selector;
      public:
        InvalidParent(Selector_Ptr parent, Backtraces traces, Selector_Ptr selector);
        virtual ~InvalidParent() throw() {};
    };

    class MissingArgument : public Base {
      protected:
        std::string fn;
        std::string arg;
        std::string fntype;
      public:
        MissingArgument(ParserState pstate, Backtraces traces, std::string fn, std::string arg, std::string fntype);
        virtual ~MissingArgument() throw() {};
    };

    class InvalidArgumentType : public Base {
      protected:
        std::string fn;
        std::string arg;
        std::string type;
        const Value_Ptr value;
      public:
        InvalidArgumentType(ParserState pstate, Backtraces traces, std::string fn, std::string arg, std::string type, const Value_Ptr value = 0);
        virtual ~InvalidArgumentType() throw() {};
    };

    class InvalidVarKwdType : public Base {
      protected:
        std::string name;
        const Argument_Ptr arg;
      public:
        InvalidVarKwdType(ParserState pstate, Backtraces traces, std::string name, const Argument_Ptr arg = 0);
        virtual ~InvalidVarKwdType() throw() {};
    };

    class InvalidSyntax : public Base {
      public:
        InvalidSyntax(ParserState pstate, Backtraces traces, std::string msg);
        virtual ~InvalidSyntax() throw() {};
    };

    class NestingLimitError : public Base {
      public:
        NestingLimitError(ParserState pstate, Backtraces traces, std::string msg = def_nesting_limit);
        virtual ~NestingLimitError() throw() {};
    };

    class DuplicateKeyError : public Base {
      protected:
        const Map& dup;
        const Expression& org;
      public:
        DuplicateKeyError(Backtraces traces, const Map& dup, const Expression& org);
        virtual const char* errtype() const { return "Error"; }
        virtual ~DuplicateKeyError() throw() {};
    };

    class TypeMismatch : public Base {
      protected:
        const Expression& var;
        const std::string type;
      public:
        TypeMismatch(Backtraces traces, const Expression& var, const std::string type);
        virtual const char* errtype() const { return "Error"; }
        virtual ~TypeMismatch() throw() {};
    };

    class InvalidValue : public Base {
      protected:
        const Expression& val;
      public:
        InvalidValue(Backtraces traces, const Expression& val);
        virtual const char* errtype() const { return "Error"; }
        virtual ~InvalidValue() throw() {};
    };

    class StackError : public Base {
      protected:
        const AST_Node& node;
      public:
        StackError(Backtraces traces, const AST_Node& node);
        virtual const char* errtype() const { return "SystemStackError"; }
        virtual ~StackError() throw() {};
    };

    /* common virtual base class (has no pstate or trace) */
    class OperationError : public std::runtime_error {
      protected:
        std::string msg;
      public:
        OperationError(std::string msg = def_op_msg)
        : std::runtime_error(msg), msg(msg)
        {};
      public:
        virtual const char* errtype() const { return "Error"; }
        virtual const char* what() const throw() { return msg.c_str(); }
        virtual ~OperationError() throw() {};
    };

    class ZeroDivisionError : public OperationError {
      protected:
        const Expression& lhs;
        const Expression& rhs;
      public:
        ZeroDivisionError(const Expression& lhs, const Expression& rhs);
        virtual const char* errtype() const { return "ZeroDivisionError"; }
        virtual ~ZeroDivisionError() throw() {};
    };

    class IncompatibleUnits : public OperationError {
      protected:
        // const Sass::UnitType lhs;
        // const Sass::UnitType rhs;
      public:
        IncompatibleUnits(const Units& lhs, const Units& rhs);
        IncompatibleUnits(const UnitType lhs, const UnitType rhs);
        virtual ~IncompatibleUnits() throw() {};
    };

    class UndefinedOperation : public OperationError {
      protected:
        Expression_Ptr_Const lhs;
        Expression_Ptr_Const rhs;
        const Sass_OP op;
      public:
        UndefinedOperation(Expression_Ptr_Const lhs, Expression_Ptr_Const rhs, enum Sass_OP op);
        // virtual const char* errtype() const { return "Error"; }
        virtual ~UndefinedOperation() throw() {};
    };

    class InvalidNullOperation : public UndefinedOperation {
      public:
        InvalidNullOperation(Expression_Ptr_Const lhs, Expression_Ptr_Const rhs, enum Sass_OP op);
        virtual ~InvalidNullOperation() throw() {};
    };

    class AlphaChannelsNotEqual : public OperationError {
      protected:
        Expression_Ptr_Const lhs;
        Expression_Ptr_Const rhs;
        const Sass_OP op;
      public:
        AlphaChannelsNotEqual(Expression_Ptr_Const lhs, Expression_Ptr_Const rhs, enum Sass_OP op);
        // virtual const char* errtype() const { return "Error"; }
        virtual ~AlphaChannelsNotEqual() throw() {};
    };

    class SassValueError : public Base {
      public:
        SassValueError(Backtraces traces, ParserState pstate, OperationError& err);
        virtual ~SassValueError() throw() {};
    };

  }

  void warn(std::string msg, ParserState pstate);
  void warn(std::string msg, ParserState pstate, Backtrace* bt);
  void warning(std::string msg, ParserState pstate);

  void deprecated_function(std::string msg, ParserState pstate);
  void deprecated(std::string msg, std::string msg2, bool with_column, ParserState pstate);
  void deprecated_bind(std::string msg, ParserState pstate);
  // void deprecated(std::string msg, ParserState pstate, Backtrace* bt);

  void coreError(std::string msg, ParserState pstate);
  void error(std::string msg, ParserState pstate, Backtraces& traces);

}

#endif
#ifndef SASS_FN_UTILS_H
#define SASS_FN_UTILS_H

namespace Sass {

  #define FN_PROTOTYPE \
    Env& env, \
    Env& d_env, \
    Context& ctx, \
    Signature sig, \
    ParserState pstate, \
    Backtraces& traces, \
    SelectorStack& selector_stack

  typedef const char* Signature;
  typedef PreValue_Ptr (*Native_Function)(FN_PROTOTYPE);
  #define BUILT_IN(name) PreValue_Ptr name(FN_PROTOTYPE)

  #define ARG(argname, argtype) get_arg<argtype>(argname, env, sig, pstate, traces)
  // special function for weird hsla percent (10px == 10% == 10 != 0.1)
  #define ARGVAL(argname) get_arg_val(argname, env, sig, pstate, traces) // double

  Definition_Ptr make_native_function(Signature, Native_Function, Context& ctx);
  Definition_Ptr make_c_function(Sass_Function_Entry c_func, Context& ctx);

  namespace Functions {

    template <typename T>
    T* get_arg(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces)
    {
      T* val = Cast<T>(env[argname]);
      if (!val) {
        error("argument `" + argname + "` of `" + sig + "` must be a " + T::type_name(), pstate, traces);
      }
      return val;
    }

    Map_Ptr get_arg_m(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces); // maps only
    Number_Ptr get_arg_n(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces); // numbers only
    double alpha_num(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces); // colors only
    double color_num(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces); // colors only
    double get_arg_r(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces, double lo, double hi); // colors only
    double get_arg_val(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces); // shared
    Selector_List_Obj get_arg_sels(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces, Context& ctx); // selectors only
    Compound_Selector_Obj get_arg_sel(const std::string& argname, Env& env, Signature sig, ParserState pstate, Backtraces traces, Context& ctx); // selectors only

  }

}

#endif

#ifndef SASS_H
#define SASS_H

// #define DEBUG 1

// include API headers
#ifndef SASS_VERSION_H
#define SASS_VERSION_H

#ifndef LIBSASS_VERSION
#define LIBSASS_VERSION "[NA]"
#endif

#ifndef LIBSASS_LANGUAGE_VERSION
#define LIBSASS_LANGUAGE_VERSION "3.5"
#endif

#endif
/**
 * sass2scss
 * Licensed under the MIT License
 * Copyright (c) Marcel Greter
 */

#ifndef SASS2SCSS_H
#define SASS2SCSS_H

#ifdef _WIN32

  /* You should define ADD_EXPORTS *only* when building the DLL. */
  #ifdef ADD_EXPORTS
    #define ADDAPI __declspec(dllexport)
	#define ADDCALL __cdecl
  #else
    #define ADDAPI
	#define ADDCALL
  #endif

#else /* _WIN32 not defined. */

  /* Define with no value on non-Windows OSes. */
  #define ADDAPI
  #define ADDCALL

#endif

#ifdef __cplusplus

#include <stack>

#ifndef SASS2SCSS_VERSION
// Hardcode once the file is copied from
// https://github.com/mgreter/sass2scss
#define SASS2SCSS_VERSION "1.1.1"
#endif

// add namespace for c++
namespace Sass
{

	// pretty print options
	const int SASS2SCSS_PRETTIFY_0 = 0;
	const int SASS2SCSS_PRETTIFY_1 = 1;
	const int SASS2SCSS_PRETTIFY_2 = 2;
	const int SASS2SCSS_PRETTIFY_3 = 3;

	// remove one-line comment
	const int SASS2SCSS_KEEP_COMMENT    =  32;
	// remove multi-line comments
	const int SASS2SCSS_STRIP_COMMENT   =  64;
	// convert one-line to multi-line
	const int SASS2SCSS_CONVERT_COMMENT = 128;

	// String for finding something interesting
	const std::string SASS2SCSS_FIND_WHITESPACE = " \t\n\v\f\r";

	// converter struct
	// holding all states
	struct converter
	{
		// bit options
		int options;
		// is selector
		bool selector;
		// concat lists
		bool comma;
		// has property
		bool property;
		// has semicolon
		bool semicolon;
		// comment context
		std::string comment;
		// flag end of file
		bool end_of_file;
		// whitespace buffer
		std::string whitespace;
		// context/block stack
		std::stack<std::string> indents;
	};

	// function only available in c++ code
	char* sass2scss (const std::string& sass, const int options);

}
// EO namespace

// declare for c
extern "C" {
#endif

	// prettyfy print options
	#define SASS2SCSS_PRETTIFY_0   0
	#define SASS2SCSS_PRETTIFY_1   1
	#define SASS2SCSS_PRETTIFY_2   2
	#define SASS2SCSS_PRETTIFY_3   3

	// keep one-line comments
	#define SASS2SCSS_KEEP_COMMENT     32
	// remove multi-line comments
	#define SASS2SCSS_STRIP_COMMENT    64
	// convert one-line to multi-line
	#define SASS2SCSS_CONVERT_COMMENT  128

	// available to c and c++ code
	ADDAPI char* ADDCALL sass2scss (const char* sass, const int options);

	// Get compiled sass2scss version
	ADDAPI const char* ADDCALL sass2scss_version(void);

#ifdef __cplusplus
} // __cplusplus defined.
#endif

#endif

#endif


namespace Sass {

  // easier to search with name
  const bool DELAYED = true;

  // ToDo: should this really be hardcoded
  // Note: most methods follow precision option
  const double NUMBER_EPSILON = 1e-12;

  // macro to test if numbers are equal within a small error margin
  #define NEAR_EQUAL(lhs, rhs) std::fabs(lhs - rhs) < NUMBER_EPSILON

  // ToDo: where does this fit best?
  // We don't share this with C-API?
  class Operand {
    public:
      Operand(Sass_OP operand, bool ws_before = false, bool ws_after = false)
      : operand(operand), ws_before(ws_before), ws_after(ws_after)
      { }
    public:
      enum Sass_OP operand;
      bool ws_before;
      bool ws_after;
  };

  //////////////////////////////////////////////////////////
  // `hash_combine` comes from boost (functional/hash):
  // http://www.boost.org/doc/libs/1_35_0/doc/html/hash/combine.html
  // Boost Software License - Version 1.0
  // http://www.boost.org/users/license.html
  template <typename T>
  void hash_combine (std::size_t& seed, const T& val)
  {
    seed ^= std::hash<T>()(val) + 0x9e3779b9
             + (seed<<6) + (seed>>2);
  }
  //////////////////////////////////////////////////////////

  const char* sass_op_to_name(enum Sass_OP op);

  const char* sass_op_separator(enum Sass_OP op);

  //////////////////////////////////////////////////////////
  // Abstract base class for all abstract syntax tree nodes.
  //////////////////////////////////////////////////////////
  class AST_Node : public SharedObj {
    ADD_PROPERTY(ParserState, pstate)
  public:
    AST_Node(ParserState pstate)
    : pstate_(pstate)
    { }
    AST_Node(const AST_Node* ptr)
    : pstate_(ptr->pstate_)
    { }

    // allow implicit conversion to string
    // needed for by SharedPtr implementation
    operator std::string() {
      return to_string();
    }

    // AST_Node(AST_Node& ptr) = delete;

    virtual ~AST_Node() = 0;
    virtual size_t hash() const { return 0; }
    virtual std::string inspect() const { return to_string({ INSPECT, 5 }); }
    virtual std::string to_sass() const { return to_string({ TO_SASS, 5 }); }
    virtual const std::string to_string(Sass_Inspect_Options opt) const;
    virtual const std::string to_string() const;
    virtual void cloneChildren() {};
    // generic find function (not fully implemented yet)
    // ToDo: add specific implementions to all children
    virtual bool find ( bool (*f)(AST_Node_Obj) ) { return f(this); };
    void update_pstate(const ParserState& pstate);
    Offset off() { return pstate(); }
    Position pos() { return pstate(); }
    ATTACH_ABSTRACT_AST_OPERATIONS(AST_Node);
    ATTACH_ABSTRACT_CRTP_PERFORM_METHODS()
  };
  inline AST_Node::~AST_Node() { }

  //////////////////////////////////////////////////////////////////////
  // define cast template now (need complete type)
  //////////////////////////////////////////////////////////////////////

  template<class T>
  T* Cast(AST_Node* ptr) {
    return ptr && typeid(T) == typeid(*ptr) ?
           static_cast<T*>(ptr) : NULL;
  };

  template<class T>
  const T* Cast(const AST_Node* ptr) {
    return ptr && typeid(T) == typeid(*ptr) ?
           static_cast<const T*>(ptr) : NULL;
  };

  //////////////////////////////////////////////////////////////////////
  // Abstract base class for expressions. This side of the AST hierarchy
  // represents elements in value contexts, which exist primarily to be
  // evaluated and returned.
  //////////////////////////////////////////////////////////////////////
  class Expression : public AST_Node {
  public:
    enum Type {
      NONE,
      BOOLEAN,
      NUMBER,
      COLOR,
      STRING,
      LIST,
      MAP,
      SELECTOR,
      NULL_VAL,
      FUNCTION_VAL,
      C_WARNING,
      C_ERROR,
      FUNCTION,
      VARIABLE,
      PARENT,
      NUM_TYPES
    };
  private:
    // expressions in some contexts shouldn't be evaluated
    ADD_PROPERTY(bool, is_delayed)
    ADD_PROPERTY(bool, is_expanded)
    ADD_PROPERTY(bool, is_interpolant)
    ADD_PROPERTY(Type, concrete_type)
  public:
    Expression(ParserState pstate, bool d = false, bool e = false, bool i = false, Type ct = NONE);
    virtual operator bool() { return true; }
    virtual ~Expression() { }
    virtual bool is_invisible() const { return false; }

    virtual std::string type() const { return ""; }
    static std::string type_name() { return ""; }

    virtual bool is_false() { return false; }
    // virtual bool is_true() { return !is_false(); }
    virtual bool operator< (const Expression& rhs) const { return false; }
    virtual bool operator== (const Expression& rhs) const { return false; }
    inline bool operator>(const Expression& rhs) const { return rhs < *this; }
    inline bool operator!=(const Expression& rhs) const { return !(rhs == *this); }
    virtual bool eq(const Expression& rhs) const { return *this == rhs; };
    virtual void set_delayed(bool delayed) { is_delayed(delayed); }
    virtual bool has_interpolant() const { return is_interpolant(); }
    virtual bool is_left_interpolant() const { return is_interpolant(); }
    virtual bool is_right_interpolant() const { return is_interpolant(); }
    ATTACH_VIRTUAL_AST_OPERATIONS(Expression);
    size_t hash() const override { return 0; }
  };

}

/////////////////////////////////////////////////////////////////////////////////////
// Hash method specializations for std::unordered_map to work with Sass::Expression
/////////////////////////////////////////////////////////////////////////////////////

namespace std {
  template<>
  struct hash<Sass::Expression_Obj>
  {
    size_t operator()(Sass::Expression_Obj s) const
    {
      return s->hash();
    }
  };
  template<>
  struct equal_to<Sass::Expression_Obj>
  {
    bool operator()( Sass::Expression_Obj lhs,  Sass::Expression_Obj rhs) const
    {
      return lhs->hash() == rhs->hash();
    }
  };
}

namespace Sass {

  /////////////////////////////////////////////////////////////////////////////
  // Mixin class for AST nodes that should behave like vectors. Uses the
  // "Template Method" design pattern to allow subclasses to adjust their flags
  // when certain objects are pushed.
  /////////////////////////////////////////////////////////////////////////////
  template <typename T>
  class Vectorized {
    std::vector<T> elements_;
  protected:
    mutable size_t hash_;
    void reset_hash() { hash_ = 0; }
    virtual void adjust_after_pushing(T element) { }
  public:
    Vectorized(size_t s = 0) : hash_(0)
    { elements_.reserve(s); }
    virtual ~Vectorized() = 0;
    size_t length() const   { return elements_.size(); }
    bool empty() const      { return elements_.empty(); }
    void clear()            { return elements_.clear(); }
    T last() const          { return elements_.back(); }
    T first() const         { return elements_.front(); }
    T& operator[](size_t i) { return elements_[i]; }
    virtual const T& at(size_t i) const { return elements_.at(i); }
    virtual T& at(size_t i) { return elements_.at(i); }
    const T& get(size_t i) const { return elements_[i]; }
    const T& operator[](size_t i) const { return elements_[i]; }
    virtual void append(T element)
    {
      if (element) {
        reset_hash();
        elements_.push_back(element);
        adjust_after_pushing(element);
      }
    }
    virtual void concat(Vectorized* v)
    {
      for (size_t i = 0, L = v->length(); i < L; ++i) this->append((*v)[i]);
    }
    Vectorized& unshift(T element)
    {
      elements_.insert(elements_.begin(), element);
      return *this;
    }
    std::vector<T>& elements() { return elements_; }
    const std::vector<T>& elements() const { return elements_; }
    std::vector<T>& elements(std::vector<T>& e) { elements_ = e; return elements_; }

    virtual size_t hash() const
    {
      if (hash_ == 0) {
        for (const T& el : elements_) {
          hash_combine(hash_, el->hash());
        }
      }
      return hash_;
    }

    template <typename P, typename V>
    typename std::vector<T>::iterator insert(P position, const V& val) {
      reset_hash();
      return elements_.insert(position, val);
    }

    typename std::vector<T>::iterator end() { return elements_.end(); }
    typename std::vector<T>::iterator begin() { return elements_.begin(); }
    typename std::vector<T>::const_iterator end() const { return elements_.end(); }
    typename std::vector<T>::const_iterator begin() const { return elements_.begin(); }
    typename std::vector<T>::iterator erase(typename std::vector<T>::iterator el) { return elements_.erase(el); }
    typename std::vector<T>::const_iterator erase(typename std::vector<T>::const_iterator el) { return elements_.erase(el); }

  };
  template <typename T>
  inline Vectorized<T>::~Vectorized() { }

  /////////////////////////////////////////////////////////////////////////////
  // Mixin class for AST nodes that should behave like a hash table. Uses an
  // extra <std::vector> internally to maintain insertion order for interation.
  /////////////////////////////////////////////////////////////////////////////
  class Hashed {
  private:
    ExpressionMap elements_;
    std::vector<Expression_Obj> list_;
  protected:
    mutable size_t hash_;
    Expression_Obj duplicate_key_;
    void reset_hash() { hash_ = 0; }
    void reset_duplicate_key() { duplicate_key_ = {}; }
    virtual void adjust_after_pushing(std::pair<Expression_Obj, Expression_Obj> p) { }
  public:
    Hashed(size_t s = 0)
    : elements_(ExpressionMap(s)),
      list_(std::vector<Expression_Obj>()),
      hash_(0), duplicate_key_({})
    { elements_.reserve(s); list_.reserve(s); }
    virtual ~Hashed();
    size_t length() const                  { return list_.size(); }
    bool empty() const                     { return list_.empty(); }
    bool has(Expression_Obj k) const          { return elements_.count(k) == 1; }
    Expression_Obj at(Expression_Obj k) const;
    bool has_duplicate_key() const         { return duplicate_key_ != 0; }
    Expression_Obj get_duplicate_key() const  { return duplicate_key_; }
    const ExpressionMap elements() { return elements_; }
    Hashed& operator<<(std::pair<Expression_Obj, Expression_Obj> p)
    {
      reset_hash();

      if (!has(p.first)) list_.push_back(p.first);
      else if (!duplicate_key_) duplicate_key_ = p.first;

      elements_[p.first] = p.second;

      adjust_after_pushing(p);
      return *this;
    }
    Hashed& operator+=(Hashed* h)
    {
      if (length() == 0) {
        this->elements_ = h->elements_;
        this->list_ = h->list_;
        return *this;
      }

      for (auto key : h->keys()) {
        *this << std::make_pair(key, h->at(key));
      }

      reset_duplicate_key();
      return *this;
    }
    const ExpressionMap& pairs() const { return elements_; }
    const std::vector<Expression_Obj>& keys() const { return list_; }

//    std::unordered_map<Expression_Obj, Expression_Obj>::iterator end() { return elements_.end(); }
//    std::unordered_map<Expression_Obj, Expression_Obj>::iterator begin() { return elements_.begin(); }
//    std::unordered_map<Expression_Obj, Expression_Obj>::const_iterator end() const { return elements_.end(); }
//    std::unordered_map<Expression_Obj, Expression_Obj>::const_iterator begin() const { return elements_.begin(); }

  };
  inline Hashed::~Hashed() { }


  /////////////////////////////////////////////////////////////////////////
  // Abstract base class for statements. This side of the AST hierarchy
  // represents elements in expansion contexts, which exist primarily to be
  // rewritten and macro-expanded.
  /////////////////////////////////////////////////////////////////////////
  class Statement : public AST_Node {
  public:
    enum Type {
      NONE,
      RULESET,
      MEDIA,
      DIRECTIVE,
      SUPPORTS,
      ATROOT,
      BUBBLE,
      CONTENT,
      KEYFRAMERULE,
      DECLARATION,
      ASSIGNMENT,
      IMPORT_STUB,
      IMPORT,
      COMMENT,
      WARNING,
      RETURN,
      EXTEND,
      ERROR,
      DEBUGSTMT,
      WHILE,
      EACH,
      FOR,
      IF
    };
  private:
    ADD_PROPERTY(Type, statement_type)
    ADD_PROPERTY(size_t, tabs)
    ADD_PROPERTY(bool, group_end)
  public:
    Statement(ParserState pstate, Type st = NONE, size_t t = 0);
    virtual ~Statement() = 0; // virtual destructor
    // needed for rearranging nested rulesets during CSS emission
    virtual bool bubbles();
    virtual bool has_content();
    virtual bool is_invisible() const;
    ATTACH_VIRTUAL_AST_OPERATIONS(Statement)
  };
  inline Statement::~Statement() { }

  ////////////////////////
  // Blocks of statements.
  ////////////////////////
  class Block final : public Statement, public Vectorized<Statement_Obj> {
    ADD_PROPERTY(bool, is_root)
    // needed for properly formatted CSS emission
  protected:
    void adjust_after_pushing(Statement_Obj s) override {}
  public:
    Block(ParserState pstate, size_t s = 0, bool r = false);
    bool has_content() override;
    ATTACH_AST_OPERATIONS(Block)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////
  // Abstract base class for statements that contain blocks of statements.
  ////////////////////////////////////////////////////////////////////////
  class Has_Block : public Statement {
    ADD_PROPERTY(Block_Obj, block)
  public:
    Has_Block(ParserState pstate, Block_Obj b);
    Has_Block(const Has_Block* ptr); // copy constructor
    virtual ~Has_Block() = 0; // virtual destructor
    virtual bool has_content() override;
  };
  inline Has_Block::~Has_Block() { }

  /////////////////////////////////////////////////////////////////////////////
  // Rulesets (i.e., sets of styles headed by a selector and containing a block
  // of style declarations.
  /////////////////////////////////////////////////////////////////////////////
  class Ruleset final : public Has_Block {
    ADD_PROPERTY(Selector_List_Obj, selector)
    ADD_PROPERTY(bool, is_root);
  public:
    Ruleset(ParserState pstate, Selector_List_Obj s = {}, Block_Obj b = {});
    bool is_invisible() const override;
    ATTACH_AST_OPERATIONS(Ruleset)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////
  // Bubble.
  /////////////////
  class Bubble final : public Statement {
    ADD_PROPERTY(Statement_Obj, node)
    ADD_PROPERTY(bool, group_end)
  public:
    Bubble(ParserState pstate, Statement_Obj n, Statement_Obj g = {}, size_t t = 0);
    bool bubbles() override;
    ATTACH_AST_OPERATIONS(Bubble)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////
  // Trace.
  /////////////////
  class Trace final : public Has_Block {
    ADD_CONSTREF(char, type)
    ADD_CONSTREF(std::string, name)
  public:
    Trace(ParserState pstate, std::string n, Block_Obj b = {}, char type = 'm');
    ATTACH_AST_OPERATIONS(Trace)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////
  // Media queries.
  /////////////////
  class Media_Block final : public Has_Block {
    ADD_PROPERTY(List_Obj, media_queries)
  public:
    Media_Block(ParserState pstate, List_Obj mqs, Block_Obj b);
    bool bubbles() override;
    bool is_invisible() const override;
    ATTACH_AST_OPERATIONS(Media_Block)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////////////////////////////////
  // At-rules -- arbitrary directives beginning with "@" that may have an
  // optional statement block.
  ///////////////////////////////////////////////////////////////////////
  class Directive final : public Has_Block {
    ADD_CONSTREF(std::string, keyword)
    ADD_PROPERTY(Selector_List_Obj, selector)
    ADD_PROPERTY(Expression_Obj, value)
  public:
    Directive(ParserState pstate, std::string kwd, Selector_List_Obj sel = {}, Block_Obj b = {}, Expression_Obj val = {});
    bool bubbles() override;
    bool is_media();
    bool is_keyframes();
    ATTACH_AST_OPERATIONS(Directive)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////////////////////////////////
  // Keyframe-rules -- the child blocks of "@keyframes" nodes.
  ///////////////////////////////////////////////////////////////////////
  class Keyframe_Rule final : public Has_Block {
    // according to css spec, this should be <keyframes-name>
    // <keyframes-name> = <custom-ident> | <string>
    ADD_PROPERTY(Selector_List_Obj, name)
  public:
    Keyframe_Rule(ParserState pstate, Block_Obj b);
    ATTACH_AST_OPERATIONS(Keyframe_Rule)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////
  // Declarations -- style rules consisting of a property name and values.
  ////////////////////////////////////////////////////////////////////////
  class Declaration final : public Has_Block {
    ADD_PROPERTY(String_Obj, property)
    ADD_PROPERTY(Expression_Obj, value)
    ADD_PROPERTY(bool, is_important)
    ADD_PROPERTY(bool, is_custom_property)
    ADD_PROPERTY(bool, is_indented)
  public:
    Declaration(ParserState pstate, String_Obj prop, Expression_Obj val, bool i = false, bool c = false, Block_Obj b = {});
    bool is_invisible() const override;
    ATTACH_AST_OPERATIONS(Declaration)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////
  // Assignments -- variable and value.
  /////////////////////////////////////
  class Assignment final : public Statement {
    ADD_CONSTREF(std::string, variable)
    ADD_PROPERTY(Expression_Obj, value)
    ADD_PROPERTY(bool, is_default)
    ADD_PROPERTY(bool, is_global)
  public:
    Assignment(ParserState pstate, std::string var, Expression_Obj val, bool is_default = false, bool is_global = false);
    ATTACH_AST_OPERATIONS(Assignment)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////////
  // Import directives. CSS and Sass import lists can be intermingled, so it's
  // necessary to store a list of each in an Import node.
  ////////////////////////////////////////////////////////////////////////////
  class Import final : public Statement {
    std::vector<Expression_Obj> urls_;
    std::vector<Include>        incs_;
    ADD_PROPERTY(List_Obj,      import_queries);
  public:
    Import(ParserState pstate);
    std::vector<Include>& incs();
    std::vector<Expression_Obj>& urls();
    ATTACH_AST_OPERATIONS(Import)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  // not yet resolved single import
  // so far we only know requested name
  class Import_Stub final : public Statement {
    Include resource_;
  public:
    Import_Stub(ParserState pstate, Include res);
    Include resource();
    std::string imp_path();
    std::string abs_path();
    ATTACH_AST_OPERATIONS(Import_Stub)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////
  // The Sass `@warn` directive.
  //////////////////////////////
  class Warning final : public Statement {
    ADD_PROPERTY(Expression_Obj, message)
  public:
    Warning(ParserState pstate, Expression_Obj msg);
    ATTACH_AST_OPERATIONS(Warning)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////
  // The Sass `@error` directive.
  ///////////////////////////////
  class Error final : public Statement {
    ADD_PROPERTY(Expression_Obj, message)
  public:
    Error(ParserState pstate, Expression_Obj msg);
    ATTACH_AST_OPERATIONS(Error)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////
  // The Sass `@debug` directive.
  ///////////////////////////////
  class Debug final : public Statement {
    ADD_PROPERTY(Expression_Obj, value)
  public:
    Debug(ParserState pstate, Expression_Obj val);
    ATTACH_AST_OPERATIONS(Debug)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////
  // CSS comments. These may be interpolated.
  ///////////////////////////////////////////
  class Comment final : public Statement {
    ADD_PROPERTY(String_Obj, text)
    ADD_PROPERTY(bool, is_important)
  public:
    Comment(ParserState pstate, String_Obj txt, bool is_important);
    virtual bool is_invisible() const override;
    ATTACH_AST_OPERATIONS(Comment)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////
  // The Sass `@if` control directive.
  ////////////////////////////////////
  class If final : public Has_Block {
    ADD_PROPERTY(Expression_Obj, predicate)
    ADD_PROPERTY(Block_Obj, alternative)
  public:
    If(ParserState pstate, Expression_Obj pred, Block_Obj con, Block_Obj alt = {});
    virtual bool has_content() override;
    ATTACH_AST_OPERATIONS(If)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////
  // The Sass `@for` control directive.
  /////////////////////////////////////
  class For final : public Has_Block {
    ADD_CONSTREF(std::string, variable)
    ADD_PROPERTY(Expression_Obj, lower_bound)
    ADD_PROPERTY(Expression_Obj, upper_bound)
    ADD_PROPERTY(bool, is_inclusive)
  public:
    For(ParserState pstate, std::string var, Expression_Obj lo, Expression_Obj hi, Block_Obj b, bool inc);
    ATTACH_AST_OPERATIONS(For)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////////////
  // The Sass `@each` control directive.
  //////////////////////////////////////
  class Each final : public Has_Block {
    ADD_PROPERTY(std::vector<std::string>, variables)
    ADD_PROPERTY(Expression_Obj, list)
  public:
    Each(ParserState pstate, std::vector<std::string> vars, Expression_Obj lst, Block_Obj b);
    ATTACH_AST_OPERATIONS(Each)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////
  // The Sass `@while` control directive.
  ///////////////////////////////////////
  class While final : public Has_Block {
    ADD_PROPERTY(Expression_Obj, predicate)
  public:
    While(ParserState pstate, Expression_Obj pred, Block_Obj b);
    ATTACH_AST_OPERATIONS(While)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////////////////
  // The @return directive for use inside SassScript functions.
  /////////////////////////////////////////////////////////////
  class Return final : public Statement {
    ADD_PROPERTY(Expression_Obj, value)
  public:
    Return(ParserState pstate, Expression_Obj val);
    ATTACH_AST_OPERATIONS(Return)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////
  // The Sass `@extend` directive.
  ////////////////////////////////
  class Extension final : public Statement {
    ADD_PROPERTY(Selector_List_Obj, selector)
  public:
    Extension(ParserState pstate, Selector_List_Obj s);
    ATTACH_AST_OPERATIONS(Extension)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////////////////////////////////
  // Definitions for both mixins and functions. The two cases are distinguished
  // by a type tag.
  /////////////////////////////////////////////////////////////////////////////
  class Definition final : public Has_Block {
  public:
    enum Type { MIXIN, FUNCTION };
    ADD_CONSTREF(std::string, name)
    ADD_PROPERTY(Parameters_Obj, parameters)
    ADD_PROPERTY(Env*, environment)
    ADD_PROPERTY(Type, type)
    ADD_PROPERTY(Native_Function, native_function)
    ADD_PROPERTY(Sass_Function_Entry, c_function)
    ADD_PROPERTY(void*, cookie)
    ADD_PROPERTY(bool, is_overload_stub)
    ADD_PROPERTY(Signature, signature)
  public:
    Definition(ParserState pstate,
               std::string n,
               Parameters_Obj params,
               Block_Obj b,
               Type t);
    Definition(ParserState pstate,
               Signature sig,
               std::string n,
               Parameters_Obj params,
               Native_Function func_ptr,
               bool overload_stub = false);
    Definition(ParserState pstate,
               Signature sig,
               std::string n,
               Parameters_Obj params,
               Sass_Function_Entry c_func);
    ATTACH_AST_OPERATIONS(Definition)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////////////
  // Mixin calls (i.e., `@include ...`).
  //////////////////////////////////////
  class Mixin_Call final : public Has_Block {
    ADD_CONSTREF(std::string, name)
    ADD_PROPERTY(Arguments_Obj, arguments)
    ADD_PROPERTY(Parameters_Obj, block_parameters)
  public:
    Mixin_Call(ParserState pstate, std::string n, Arguments_Obj args, Parameters_Obj b_params = {}, Block_Obj b = {});
    ATTACH_AST_OPERATIONS(Mixin_Call)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////////////
  // The @content directive for mixin content blocks.
  ///////////////////////////////////////////////////
  class Content final : public Statement {
    ADD_PROPERTY(Arguments_Obj, arguments)
  public:
    Content(ParserState pstate, Arguments_Obj args);
    ATTACH_AST_OPERATIONS(Content)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////////
  // Arithmetic negation (logical negation is just an ordinary function call).
  ////////////////////////////////////////////////////////////////////////////
  class Unary_Expression final : public Expression {
  public:
    enum Type { PLUS, MINUS, NOT, SLASH };
  private:
    HASH_PROPERTY(Type, optype)
    HASH_PROPERTY(Expression_Obj, operand)
    mutable size_t hash_;
  public:
    Unary_Expression(ParserState pstate, Type t, Expression_Obj o);
    const std::string type_name();
    virtual bool operator==(const Expression& rhs) const override;
    size_t hash() const override;
    ATTACH_AST_OPERATIONS(Unary_Expression)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////
  // Individual argument objects for mixin and function calls.
  ////////////////////////////////////////////////////////////
  class Argument final : public Expression {
    HASH_PROPERTY(Expression_Obj, value)
    HASH_CONSTREF(std::string, name)
    ADD_PROPERTY(bool, is_rest_argument)
    ADD_PROPERTY(bool, is_keyword_argument)
    mutable size_t hash_;
  public:
    Argument(ParserState pstate, Expression_Obj val, std::string n = "", bool rest = false, bool keyword = false);
    void set_delayed(bool delayed) override;
    bool operator==(const Expression& rhs) const override;
    size_t hash() const override;
    ATTACH_AST_OPERATIONS(Argument)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////
  // Argument lists -- in their own class to facilitate context-sensitive
  // error checking (e.g., ensuring that all ordinal arguments precede all
  // named arguments).
  ////////////////////////////////////////////////////////////////////////
  class Arguments final : public Expression, public Vectorized<Argument_Obj> {
    ADD_PROPERTY(bool, has_named_arguments)
    ADD_PROPERTY(bool, has_rest_argument)
    ADD_PROPERTY(bool, has_keyword_argument)
  protected:
    void adjust_after_pushing(Argument_Obj a) override;
  public:
    Arguments(ParserState pstate);
    void set_delayed(bool delayed) override;
    Argument_Obj get_rest_argument();
    Argument_Obj get_keyword_argument();
    ATTACH_AST_OPERATIONS(Arguments)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////
  // Media queries.
  /////////////////
  class Media_Query final : public Expression,
                            public Vectorized<Media_Query_Expression_Obj> {
    ADD_PROPERTY(String_Obj, media_type)
    ADD_PROPERTY(bool, is_negated)
    ADD_PROPERTY(bool, is_restricted)
  public:
    Media_Query(ParserState pstate, String_Obj t = {}, size_t s = 0, bool n = false, bool r = false);
    ATTACH_AST_OPERATIONS(Media_Query)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////
  // Media expressions (for use inside media queries).
  ////////////////////////////////////////////////////
  class Media_Query_Expression final : public Expression {
    ADD_PROPERTY(Expression_Obj, feature)
    ADD_PROPERTY(Expression_Obj, value)
    ADD_PROPERTY(bool, is_interpolated)
  public:
    Media_Query_Expression(ParserState pstate, Expression_Obj f, Expression_Obj v, bool i = false);
    ATTACH_AST_OPERATIONS(Media_Query_Expression)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////
  // At root expressions (for use inside @at-root).
  /////////////////////////////////////////////////
  class At_Root_Query final : public Expression {
  private:
    ADD_PROPERTY(Expression_Obj, feature)
    ADD_PROPERTY(Expression_Obj, value)
  public:
    At_Root_Query(ParserState pstate, Expression_Obj f = {}, Expression_Obj v = {}, bool i = false);
    bool exclude(std::string str);
    ATTACH_AST_OPERATIONS(At_Root_Query)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////
  // At-root.
  ///////////
  class At_Root_Block final : public Has_Block {
    ADD_PROPERTY(At_Root_Query_Obj, expression)
  public:
    At_Root_Block(ParserState pstate, Block_Obj b = {}, At_Root_Query_Obj e = {});
    bool bubbles() override;
    bool exclude_node(Statement_Obj s);
    ATTACH_AST_OPERATIONS(At_Root_Block)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////////////
  // Individual parameter objects for mixins and functions.
  /////////////////////////////////////////////////////////
  class Parameter final : public AST_Node {
    ADD_CONSTREF(std::string, name)
    ADD_PROPERTY(Expression_Obj, default_value)
    ADD_PROPERTY(bool, is_rest_parameter)
  public:
    Parameter(ParserState pstate, std::string n, Expression_Obj def = {}, bool rest = false);
    ATTACH_AST_OPERATIONS(Parameter)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////////////////////////////
  // Parameter lists -- in their own class to facilitate context-sensitive
  // error checking (e.g., ensuring that all optional parameters follow all
  // required parameters).
  /////////////////////////////////////////////////////////////////////////
  class Parameters final : public AST_Node, public Vectorized<Parameter_Obj> {
    ADD_PROPERTY(bool, has_optional_parameters)
    ADD_PROPERTY(bool, has_rest_parameter)
  protected:
    void adjust_after_pushing(Parameter_Obj p) override;
  public:
    Parameters(ParserState pstate);
    ATTACH_AST_OPERATIONS(Parameters)
    ATTACH_CRTP_PERFORM_METHODS()
  };

}

#ifndef SASS_AST_VALUES_H
#define SASS_AST_VALUES_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.




namespace Sass {

  //////////////////////////////////////////////////////////////////////
  // Still just an expression, but with a to_string method
  //////////////////////////////////////////////////////////////////////
  class PreValue : public Expression {
  public:
    PreValue(ParserState pstate, bool d = false, bool e = false, bool i = false, Type ct = NONE);
    ATTACH_VIRTUAL_AST_OPERATIONS(PreValue);
    virtual ~PreValue() { }
  };

  //////////////////////////////////////////////////////////////////////
  // base class for values that support operations
  //////////////////////////////////////////////////////////////////////
  class Value : public PreValue {
  public:
    Value(ParserState pstate, bool d = false, bool e = false, bool i = false, Type ct = NONE);
    ATTACH_VIRTUAL_AST_OPERATIONS(Value);
    virtual bool operator== (const Expression& rhs) const override = 0;
  };

  ///////////////////////////////////////////////////////////////////////
  // Lists of values, both comma- and space-separated (distinguished by a
  // type-tag.) Also used to represent variable-length argument lists.
  ///////////////////////////////////////////////////////////////////////
  class List : public Value, public Vectorized<Expression_Obj> {
    void adjust_after_pushing(Expression_Obj e) override { is_expanded(false); }
  private:
    ADD_PROPERTY(enum Sass_Separator, separator)
    ADD_PROPERTY(bool, is_arglist)
    ADD_PROPERTY(bool, is_bracketed)
    ADD_PROPERTY(bool, from_selector)
  public:
    List(ParserState pstate, size_t size = 0, enum Sass_Separator sep = SASS_SPACE, bool argl = false, bool bracket = false);
    std::string type() const override { return is_arglist_ ? "arglist" : "list"; }
    static std::string type_name() { return "list"; }
    const char* sep_string(bool compressed = false) const {
      return separator() == SASS_SPACE ?
        " " : (compressed ? "," : ", ");
    }
    bool is_invisible() const override { return empty() && !is_bracketed(); }
    Expression_Obj value_at_index(size_t i);

    virtual size_t hash() const override;
    virtual size_t size() const;
    virtual void set_delayed(bool delayed) override;
    virtual bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(List)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////////////////////////////////
  // Key value paris.
  ///////////////////////////////////////////////////////////////////////
  class Map : public Value, public Hashed {
    void adjust_after_pushing(std::pair<Expression_Obj, Expression_Obj> p) override { is_expanded(false); }
  public:
    Map(ParserState pstate, size_t size = 0);
    std::string type() const override { return "map"; }
    static std::string type_name() { return "map"; }
    bool is_invisible() const override { return empty(); }
    List_Obj to_list(ParserState& pstate);

    virtual size_t hash() const override;
    virtual bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(Map)
    ATTACH_CRTP_PERFORM_METHODS()
  };


  //////////////////////////////////////////////////////////////////////////
  // Binary expressions. Represents logical, relational, and arithmetic
  // operations. Templatized to avoid large switch statements and repetitive
  // subclassing.
  //////////////////////////////////////////////////////////////////////////
  class Binary_Expression : public PreValue {
  private:
    HASH_PROPERTY(Operand, op)
    HASH_PROPERTY(Expression_Obj, left)
    HASH_PROPERTY(Expression_Obj, right)
    mutable size_t hash_;
  public:
    Binary_Expression(ParserState pstate,
                      Operand op, Expression_Obj lhs, Expression_Obj rhs);

    const std::string type_name();
    const std::string separator();
    bool is_left_interpolant(void) const override;
    bool is_right_interpolant(void) const override;
    bool has_interpolant() const override;

    virtual void set_delayed(bool delayed) override;

    virtual bool operator==(const Expression& rhs) const override;

    virtual size_t hash() const override;
    enum Sass_OP optype() const { return op_.operand; }
    ATTACH_AST_OPERATIONS(Binary_Expression)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////
  // Function reference.
  ////////////////////////////////////////////////////
  class Function final : public Value {
  public:
    ADD_PROPERTY(Definition_Obj, definition)
    ADD_PROPERTY(bool, is_css)
  public:
    Function(ParserState pstate, Definition_Obj def, bool css);

    std::string type() const override { return "function"; }
    static std::string type_name() { return "function"; }
    bool is_invisible() const override { return true; }

    std::string name();

    bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(Function)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////
  // Function calls.
  //////////////////
  class Function_Call final : public PreValue {
    HASH_CONSTREF(String_Obj, sname)
    HASH_PROPERTY(Arguments_Obj, arguments)
    HASH_PROPERTY(Function_Obj, func)
    ADD_PROPERTY(bool, via_call)
    ADD_PROPERTY(void*, cookie)
    mutable size_t hash_;
  public:
    Function_Call(ParserState pstate, std::string n, Arguments_Obj args, void* cookie);
    Function_Call(ParserState pstate, std::string n, Arguments_Obj args, Function_Obj func);
    Function_Call(ParserState pstate, std::string n, Arguments_Obj args);

    Function_Call(ParserState pstate, String_Obj n, Arguments_Obj args, void* cookie);
    Function_Call(ParserState pstate, String_Obj n, Arguments_Obj args, Function_Obj func);
    Function_Call(ParserState pstate, String_Obj n, Arguments_Obj args);

    std::string name() const;
    bool is_css();

    bool operator==(const Expression& rhs) const override;

    size_t hash() const override;

    ATTACH_AST_OPERATIONS(Function_Call)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////
  // Variable references.
  ///////////////////////
  class Variable final : public PreValue {
    ADD_CONSTREF(std::string, name)
  public:
    Variable(ParserState pstate, std::string n);
    virtual bool operator==(const Expression& rhs) const override;
    virtual size_t hash() const override;
    ATTACH_AST_OPERATIONS(Variable)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////
  // Numbers, percentages, dimensions, and colors.
  ////////////////////////////////////////////////
  class Number final : public Value, public Units {
    HASH_PROPERTY(double, value)
    ADD_PROPERTY(bool, zero)
    mutable size_t hash_;
  public:
    Number(ParserState pstate, double val, std::string u = "", bool zero = true);

    bool zero() { return zero_; }

    std::string type() const override { return "number"; }
    static std::string type_name() { return "number"; }

    // cancel out unnecessary units
    // result will be in input units
    void reduce();

    // normalize units to defaults
    // needed to compare two numbers
    void normalize();

    size_t hash() const override;

    bool operator< (const Number& rhs) const;
    bool operator== (const Number& rhs) const;
    bool operator== (const Expression& rhs) const override;
    ATTACH_AST_OPERATIONS(Number)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////
  // Colors.
  //////////
  class Color : public Value {
    ADD_CONSTREF(std::string, disp)
    HASH_PROPERTY(double, a)
  protected:
    mutable size_t hash_;
  public:
    Color(ParserState pstate, double a = 1, const std::string disp = "");

    std::string type() const override { return "color"; }
    static std::string type_name() { return "color"; }

    virtual size_t hash() const override = 0;

    bool operator== (const Expression& rhs) const override;

    virtual Color_RGBA_Ptr toRGBA(bool copy = false) = 0;
    virtual Color_HSLA_Ptr toHSLA(bool copy = false) = 0;

    ATTACH_VIRTUAL_AST_OPERATIONS(Color)
  };

  //////////
  // Colors.
  //////////
  class Color_RGBA final : public Color {
    HASH_PROPERTY(double, r)
    HASH_PROPERTY(double, g)
    HASH_PROPERTY(double, b)
  public:
    Color_RGBA(ParserState pstate, double r, double g, double b, double a = 1, const std::string disp = "");

    std::string type() const override { return "color"; }
    static std::string type_name() { return "color"; }

    size_t hash() const override;
    Color_RGBA_Ptr toRGBA(bool copy = false) override;
    Color_HSLA_Ptr toHSLA(bool copy = false) override;

    bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(Color_RGBA)
    ATTACH_CRTP_PERFORM_METHODS()
  };


  //////////
  // Colors.
  //////////
  class Color_HSLA final : public Color {
    HASH_PROPERTY(double, h)
    HASH_PROPERTY(double, s)
    HASH_PROPERTY(double, l)
  public:
    Color_HSLA(ParserState pstate, double h, double s, double l, double a = 1, const std::string disp = "");

    std::string type() const override { return "color"; }
    static std::string type_name() { return "color"; }

    size_t hash() const override;
    Color_RGBA_Ptr toRGBA(bool copy = false) override;
    Color_HSLA_Ptr toHSLA(bool copy = false) override;

    bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(Color_HSLA)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////
  // Errors from Sass_Values.
  //////////////////////////////
  class Custom_Error final : public Value {
    ADD_CONSTREF(std::string, message)
  public:
    Custom_Error(ParserState pstate, std::string msg);
    bool operator== (const Expression& rhs) const override;
    ATTACH_AST_OPERATIONS(Custom_Error)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////
  // Warnings from Sass_Values.
  //////////////////////////////
  class Custom_Warning final : public Value {
    ADD_CONSTREF(std::string, message)
  public:
    Custom_Warning(ParserState pstate, std::string msg);
    bool operator== (const Expression& rhs) const override;
    ATTACH_AST_OPERATIONS(Custom_Warning)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////
  // Booleans.
  ////////////
  class Boolean final : public Value {
    HASH_PROPERTY(bool, value)
    mutable size_t hash_;
  public:
    Boolean(ParserState pstate, bool val);
    operator bool() override { return value_; }

    std::string type() const override { return "bool"; }
    static std::string type_name() { return "bool"; }

    size_t hash() const override;

    bool is_false() override { return !value_; }

    bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(Boolean)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////
  // Abstract base class for Sass string values. Includes interpolated and
  // "flat" strings.
  ////////////////////////////////////////////////////////////////////////
  class String : public Value {
  public:
    String(ParserState pstate, bool delayed = false);
    static std::string type_name() { return "string"; }
    virtual ~String() = 0;
    virtual void rtrim() = 0;
    virtual bool operator<(const Expression& rhs) const override {
      return this->to_string() < rhs.to_string();
    };
    ATTACH_VIRTUAL_AST_OPERATIONS(String);
    ATTACH_CRTP_PERFORM_METHODS()
  };
  inline String::~String() { };

  ///////////////////////////////////////////////////////////////////////
  // Interpolated strings. Meant to be reduced to flat strings during the
  // evaluation phase.
  ///////////////////////////////////////////////////////////////////////
  class String_Schema final : public String, public Vectorized<PreValue_Obj> {
    ADD_PROPERTY(bool, css)
    mutable size_t hash_;
  public:
    String_Schema(ParserState pstate, size_t size = 0, bool css = true);

    std::string type() const override { return "string"; }
    static std::string type_name() { return "string"; }

    bool is_left_interpolant(void) const override;
    bool is_right_interpolant(void) const override;

    bool has_interpolants();
    void rtrim() override;
    size_t hash() const override;
    virtual void set_delayed(bool delayed) override;

    bool operator==(const Expression& rhs) const override;
    ATTACH_AST_OPERATIONS(String_Schema)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////
  // Flat strings -- the lowest level of raw textual data.
  ////////////////////////////////////////////////////////
  class String_Constant : public String {
    ADD_PROPERTY(char, quote_mark)
    ADD_PROPERTY(bool, can_compress_whitespace)
    HASH_CONSTREF(std::string, value)
  protected:
    mutable size_t hash_;
  public:
    String_Constant(ParserState pstate, std::string val, bool css = true);
    String_Constant(ParserState pstate, const char* beg, bool css = true);
    String_Constant(ParserState pstate, const char* beg, const char* end, bool css = true);
    String_Constant(ParserState pstate, const Token& tok, bool css = true);
    std::string type() const override { return "string"; }
    static std::string type_name() { return "string"; }
    bool is_invisible() const override;
    virtual void rtrim() override;
    size_t hash() const override;
    bool operator==(const Expression& rhs) const override;
    // quotes are forced on inspection
    virtual std::string inspect() const override;
    ATTACH_AST_OPERATIONS(String_Constant)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////
  // Possibly quoted string (unquote on instantiation)
  ////////////////////////////////////////////////////////
  class String_Quoted final : public String_Constant {
  public:
    String_Quoted(ParserState pstate, std::string val, char q = 0,
      bool keep_utf8_escapes = false, bool skip_unquoting = false,
      bool strict_unquoting = true, bool css = true);
    bool operator==(const Expression& rhs) const override;
    // quotes are forced on inspection
    std::string inspect() const override;
    ATTACH_AST_OPERATIONS(String_Quoted)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////
  // The null value.
  //////////////////
  class Null final : public Value {
  public:
    Null(ParserState pstate);
    std::string type() const override { return "null"; }
    static std::string type_name() { return "null"; }
    bool is_invisible() const override { return true; }
    operator bool() override { return false; }
    bool is_false() override { return true; }

    size_t hash() const override;

    bool operator== (const Expression& rhs) const override;

    ATTACH_AST_OPERATIONS(Null)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////////
  // The Parent Reference Expression.
  //////////////////////////////////
  class Parent_Reference final : public Value {
  public:
    Parent_Reference(ParserState pstate);
    std::string type() const override { return "parent"; }
    static std::string type_name() { return "parent"; }
    bool operator==(const Expression& rhs) const override {
      return true; // can they ever be not equal?
    };
    ATTACH_AST_OPERATIONS(Parent_Reference)
    ATTACH_CRTP_PERFORM_METHODS()
  };

}

#endif
#ifndef SASS_AST_SUPPORTS_H
#define SASS_AST_SUPPORTS_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.




namespace Sass {

  ////////////////////
  // `@supports` rule.
  ////////////////////
  class Supports_Block : public Has_Block {
    ADD_PROPERTY(Supports_Condition_Obj, condition)
  public:
    Supports_Block(ParserState pstate, Supports_Condition_Obj condition, Block_Obj block = {});
    bool bubbles() override;
    ATTACH_AST_OPERATIONS(Supports_Block)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////////////////////////////
  // The abstract superclass of all Supports conditions.
  //////////////////////////////////////////////////////
  class Supports_Condition : public Expression {
  public:
    Supports_Condition(ParserState pstate);
    virtual bool needs_parens(Supports_Condition_Obj cond) const { return false; }
    ATTACH_AST_OPERATIONS(Supports_Condition)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////
  // An operator condition (e.g. `CONDITION1 and CONDITION2`).
  ////////////////////////////////////////////////////////////
  class Supports_Operator : public Supports_Condition {
  public:
    enum Operand { AND, OR };
  private:
    ADD_PROPERTY(Supports_Condition_Obj, left);
    ADD_PROPERTY(Supports_Condition_Obj, right);
    ADD_PROPERTY(Operand, operand);
  public:
    Supports_Operator(ParserState pstate, Supports_Condition_Obj l, Supports_Condition_Obj r, Operand o);
    virtual bool needs_parens(Supports_Condition_Obj cond) const override;
    ATTACH_AST_OPERATIONS(Supports_Operator)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////////////////
  // A negation condition (`not CONDITION`).
  //////////////////////////////////////////
  class Supports_Negation : public Supports_Condition {
  private:
    ADD_PROPERTY(Supports_Condition_Obj, condition);
  public:
    Supports_Negation(ParserState pstate, Supports_Condition_Obj c);
    virtual bool needs_parens(Supports_Condition_Obj cond) const override;
    ATTACH_AST_OPERATIONS(Supports_Negation)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////////
  // A declaration condition (e.g. `(feature: value)`).
  /////////////////////////////////////////////////////
  class Supports_Declaration : public Supports_Condition {
  private:
    ADD_PROPERTY(Expression_Obj, feature);
    ADD_PROPERTY(Expression_Obj, value);
  public:
    Supports_Declaration(ParserState pstate, Expression_Obj f, Expression_Obj v);
    virtual bool needs_parens(Supports_Condition_Obj cond) const override;
    ATTACH_AST_OPERATIONS(Supports_Declaration)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////////
  // An interpolation condition (e.g. `#{$var}`).
  ///////////////////////////////////////////////
  class Supports_Interpolation : public Supports_Condition {
  private:
    ADD_PROPERTY(Expression_Obj, value);
  public:
    Supports_Interpolation(ParserState pstate, Expression_Obj v);
    virtual bool needs_parens(Supports_Condition_Obj cond) const override;
    ATTACH_AST_OPERATIONS(Supports_Interpolation)
    ATTACH_CRTP_PERFORM_METHODS()
  };

}

#endif
#ifndef SASS_AST_SEL_H
#define SASS_AST_SEL_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.




namespace Sass {

  /////////////////////////////////////////
  // Abstract base class for CSS selectors.
  /////////////////////////////////////////
  class Selector : public Expression {
    // line break before list separator
    ADD_PROPERTY(bool, has_line_feed)
    // line break after list separator
    ADD_PROPERTY(bool, has_line_break)
    // maybe we have optional flag
    ADD_PROPERTY(bool, is_optional)
    // must not be a reference counted object
    // otherwise we create circular references
    ADD_PROPERTY(Media_Block_Ptr, media_block)
  protected:
    mutable size_t hash_;
  public:
    Selector(ParserState pstate);
    virtual ~Selector() = 0;
    size_t hash() const override = 0;
    virtual unsigned long specificity() const = 0;
    virtual int unification_order() const = 0;
    virtual void set_media_block(Media_Block_Ptr mb);
    virtual bool has_parent_ref() const;
    virtual bool has_real_parent_ref() const;
    // dispatch to correct handlers
    virtual bool operator<(const Selector& rhs) const = 0;
    virtual bool operator==(const Selector& rhs) const = 0;
    bool operator>(const Selector& rhs) const { return rhs < *this; };
    bool operator!=(const Selector& rhs) const { return !(rhs == *this); };
    ATTACH_VIRTUAL_AST_OPERATIONS(Selector);
  };
  inline Selector::~Selector() { }

  /////////////////////////////////////////////////////////////////////////
  // Interpolated selectors -- the interpolated String will be expanded and
  // re-parsed into a normal selector class.
  /////////////////////////////////////////////////////////////////////////
  class Selector_Schema final : public AST_Node {
    ADD_PROPERTY(String_Obj, contents)
    ADD_PROPERTY(bool, connect_parent);
    // must not be a reference counted object
    // otherwise we create circular references
    ADD_PROPERTY(Media_Block_Ptr, media_block)
    // store computed hash
    mutable size_t hash_;
  public:
    Selector_Schema(ParserState pstate, String_Obj c);
    bool has_parent_ref() const;
    bool has_real_parent_ref() const;
    bool operator<(const Selector& rhs) const;
    bool operator==(const Selector& rhs) const;
    // selector schema is not yet a final selector, so we do not
    // have a specificity for it yet. We need to
    virtual unsigned long specificity() const;
    size_t hash() const override;
    ATTACH_AST_OPERATIONS(Selector_Schema)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////
  // Abstract base class for simple selectors.
  ////////////////////////////////////////////
  class Simple_Selector : public Selector {
  public:
    enum Simple_Type {
      ID_SEL,
      TYPE_SEL,
      CLASS_SEL,
      PSEUDO_SEL,
      PARENT_SEL,
      WRAPPED_SEL,
      ATTRIBUTE_SEL,
      PLACEHOLDER_SEL,
    };
  public:
    HASH_CONSTREF(std::string, ns)
    HASH_CONSTREF(std::string, name)
    ADD_PROPERTY(Simple_Type, simple_type)
    HASH_PROPERTY(bool, has_ns)
  public:
    Simple_Selector(ParserState pstate, std::string n = "");
    virtual std::string ns_name() const;
    size_t hash() const override;
    bool empty() const;
    // namespace compare functions
    bool is_ns_eq(const Simple_Selector& r) const;
    // namespace query functions
    bool is_universal_ns() const;
    bool is_empty_ns() const;
    bool has_empty_ns() const;
    bool has_qualified_ns() const;
    // name query functions
    bool is_universal() const;
    virtual bool has_placeholder();

    virtual ~Simple_Selector() = 0;
    virtual Compound_Selector_Ptr unify_with(Compound_Selector_Ptr);

    virtual bool has_parent_ref() const override;
    virtual bool has_real_parent_ref() const override;
    virtual bool is_pseudo_element() const;
    virtual bool is_superselector_of(Compound_Selector_Ptr_Const sub) const;

    bool operator<(const Selector& rhs) const final override;
    bool operator==(const Selector& rhs) const final override;
    virtual bool operator<(const Selector_List& rhs) const;
    virtual bool operator==(const Selector_List& rhs) const;
    virtual bool operator<(const Complex_Selector& rhs) const;
    virtual bool operator==(const Complex_Selector& rhs) const;
    virtual bool operator<(const Compound_Selector& rhs) const;
    virtual bool operator==(const Compound_Selector& rhs) const;
    virtual bool operator<(const Simple_Selector& rhs) const;
    virtual bool operator==(const Simple_Selector& rhs) const;

    ATTACH_VIRTUAL_AST_OPERATIONS(Simple_Selector);
    ATTACH_CRTP_PERFORM_METHODS();

  };
  inline Simple_Selector::~Simple_Selector() { }

  //////////////////////////////////
  // The Parent Selector Expression.
  //////////////////////////////////
  class Parent_Selector final : public Simple_Selector {
    // a real parent selector is given by the user
    // others are added implicitly to connect the
    // selector scopes automatically when rendered
    // a Parent_Reference is never seen in selectors
    // and is only used in values (e.g. `prop: #{&};`)
    ADD_PROPERTY(bool, real)
  public:
    Parent_Selector(ParserState pstate, bool r = true);

    virtual bool has_parent_ref() const override;
    virtual bool has_real_parent_ref() const override;

    virtual unsigned long specificity() const override;
    int unification_order() const override
    {
      throw std::runtime_error("unification_order for Parent_Selector is undefined");
    }
    std::string type() const override { return "selector"; }
    static std::string type_name() { return "selector"; }
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Parent_Selector& rhs) const;
    bool operator==(const Parent_Selector& rhs) const;
    ATTACH_AST_OPERATIONS(Parent_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };


  /////////////////////////////////////////////////////////////////////////
  // Placeholder selectors (e.g., "%foo") for use in extend-only selectors.
  /////////////////////////////////////////////////////////////////////////
  class Placeholder_Selector final : public Simple_Selector {
  public:
    Placeholder_Selector(ParserState pstate, std::string n);

    int unification_order() const override
    {
      return Constants::UnificationOrder_Placeholder;
    }
    virtual ~Placeholder_Selector() {};
    virtual unsigned long specificity() const override;
    virtual bool has_placeholder() override;
    bool operator<(const Simple_Selector& rhs) const override;
    bool operator==(const Simple_Selector& rhs) const override;
    bool operator<(const Placeholder_Selector& rhs) const;
    bool operator==(const Placeholder_Selector& rhs) const;
    ATTACH_AST_OPERATIONS(Placeholder_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////////////////////////
  // Type selectors (and the universal selector) -- e.g., div, span, *.
  /////////////////////////////////////////////////////////////////////
  class Type_Selector final : public Simple_Selector {
  public:
    Type_Selector(ParserState pstate, std::string n);
    virtual unsigned long specificity() const override;
    int unification_order() const override
    {
      return Constants::UnificationOrder_Element;
    }
    Simple_Selector_Ptr unify_with(Simple_Selector_Ptr);
    Compound_Selector_Ptr unify_with(Compound_Selector_Ptr) override;
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Type_Selector& rhs) const;
    bool operator==(const Type_Selector& rhs) const;
    ATTACH_AST_OPERATIONS(Type_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////
  // Class selectors  -- i.e., .foo.
  ////////////////////////////////////////////////
  class Class_Selector final : public Simple_Selector {
  public:
    Class_Selector(ParserState pstate, std::string n);
    virtual unsigned long specificity() const override;
    int unification_order() const override
    {
      return Constants::UnificationOrder_Class;
    }
    Compound_Selector_Ptr unify_with(Compound_Selector_Ptr) override;
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Class_Selector& rhs) const;
    bool operator==(const Class_Selector& rhs) const;
    ATTACH_AST_OPERATIONS(Class_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////
  // ID selectors -- i.e., #foo.
  ////////////////////////////////////////////////
  class Id_Selector final : public Simple_Selector {
  public:
    Id_Selector(ParserState pstate, std::string n);
    virtual unsigned long specificity() const override;
    int unification_order() const override
    {
      return Constants::UnificationOrder_Id;
    }
    Compound_Selector_Ptr unify_with(Compound_Selector_Ptr) override;
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Id_Selector& rhs) const;
    bool operator==(const Id_Selector& rhs) const;
    ATTACH_AST_OPERATIONS(Id_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////////////////////
  // Attribute selectors -- e.g., [src*=".jpg"], etc.
  ///////////////////////////////////////////////////
  class Attribute_Selector final : public Simple_Selector {
    ADD_CONSTREF(std::string, matcher)
    // this cannot be changed to obj atm!!!!!!????!!!!!!!
    ADD_PROPERTY(String_Obj, value) // might be interpolated
    ADD_PROPERTY(char, modifier);
  public:
    Attribute_Selector(ParserState pstate, std::string n, std::string m, String_Obj v, char o = 0);
    size_t hash() const override;
    virtual unsigned long specificity() const override;
    int unification_order() const override
    {
      return Constants::UnificationOrder_Attribute;
    }
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Attribute_Selector& rhs) const;
    bool operator==(const Attribute_Selector& rhs) const;
    ATTACH_AST_OPERATIONS(Attribute_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  //////////////////////////////////////////////////////////////////
  // Pseudo selectors -- e.g., :first-child, :nth-of-type(...), etc.
  //////////////////////////////////////////////////////////////////
  /* '::' starts a pseudo-element, ':' a pseudo-class */
  /* Except :first-line, :first-letter, :before and :after */
  /* Note that pseudo-elements are restricted to one per selector */
  /* and occur only in the last simple_selector_sequence. */
  inline bool is_pseudo_class_element(const std::string& name)
  {
    return name == ":before"       ||
           name == ":after"        ||
           name == ":first-line"   ||
           name == ":first-letter";
  }

  // Pseudo Selector cannot have any namespace?
  class Pseudo_Selector final : public Simple_Selector {
    ADD_PROPERTY(String_Obj, expression)
  public:
    Pseudo_Selector(ParserState pstate, std::string n, String_Obj expr = {});
    virtual bool is_pseudo_element() const override;
    size_t hash() const override;
    virtual unsigned long specificity() const override;
    int unification_order() const override
    {
      if (is_pseudo_element())
        return Constants::UnificationOrder_PseudoElement;
      return Constants::UnificationOrder_PseudoClass;
    }
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Pseudo_Selector& rhs) const;
    bool operator==(const Pseudo_Selector& rhs) const;
    Compound_Selector_Ptr unify_with(Compound_Selector_Ptr) override;
    ATTACH_AST_OPERATIONS(Pseudo_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  /////////////////////////////////////////////////
  // Wrapped selector -- pseudo selector that takes a list of selectors as argument(s) e.g., :not(:first-of-type), :-moz-any(ol p.blah, ul, menu, dir)
  /////////////////////////////////////////////////
  class Wrapped_Selector final : public Simple_Selector {
    ADD_PROPERTY(Selector_List_Obj, selector)
  public:
    Wrapped_Selector(ParserState pstate, std::string n, Selector_List_Obj sel);
    using Simple_Selector::is_superselector_of;
    bool is_superselector_of(Wrapped_Selector_Ptr_Const sub) const;
    // Selectors inside the negation pseudo-class are counted like any
    // other, but the negation itself does not count as a pseudo-class.
    size_t hash() const override;
    bool has_parent_ref() const override;
    bool has_real_parent_ref() const override;
    unsigned long specificity() const override;
    int unification_order() const override
    {
      return Constants::UnificationOrder_Wrapped;
    }
    bool find ( bool (*f)(AST_Node_Obj) ) override;
    bool operator<(const Simple_Selector& rhs) const final override;
    bool operator==(const Simple_Selector& rhs) const final override;
    bool operator<(const Wrapped_Selector& rhs) const;
    bool operator==(const Wrapped_Selector& rhs) const;
    void cloneChildren() override;
    ATTACH_AST_OPERATIONS(Wrapped_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////////
  // Simple selector sequences. Maintains flags indicating whether it contains
  // any parent references or placeholders, to simplify expansion.
  ////////////////////////////////////////////////////////////////////////////
  class Compound_Selector final : public Selector, public Vectorized<Simple_Selector_Obj> {
  private:
    ComplexSelectorSet sources_;
    ADD_PROPERTY(bool, extended);
    ADD_PROPERTY(bool, has_parent_reference);
  protected:
    void adjust_after_pushing(Simple_Selector_Obj s) override
    {
      // if (s->has_reference())   has_reference(true);
      // if (s->has_placeholder()) has_placeholder(true);
    }
  public:
    Compound_Selector(ParserState pstate, size_t s = 0);
    bool contains_placeholder();
    void append(Simple_Selector_Obj element) override;
    bool is_universal() const;
    Complex_Selector_Obj to_complex();
    Compound_Selector_Ptr unify_with(Compound_Selector_Ptr rhs);
    // virtual Placeholder_Selector_Ptr find_placeholder();
    bool has_parent_ref() const override;
    bool has_real_parent_ref() const override;
    Simple_Selector_Ptr base() const;
    bool is_superselector_of(Compound_Selector_Ptr_Const sub, std::string wrapped = "") const;
    bool is_superselector_of(Complex_Selector_Ptr_Const sub, std::string wrapped = "") const;
    bool is_superselector_of(Selector_List_Ptr_Const sub, std::string wrapped = "") const;
    size_t hash() const override;
    virtual unsigned long specificity() const override;
    virtual bool has_placeholder();
    bool is_empty_reference();
    int unification_order() const override
    {
      throw std::runtime_error("unification_order for Compound_Selector is undefined");
    }
    bool find ( bool (*f)(AST_Node_Obj) ) override;

    bool operator<(const Selector& rhs) const override;
    bool operator==(const Selector& rhs) const override;
    bool operator<(const Selector_List& rhs) const;
    bool operator==(const Selector_List& rhs) const;
    bool operator<(const Complex_Selector& rhs) const;
    bool operator==(const Complex_Selector& rhs) const;
    bool operator<(const Compound_Selector& rhs) const;
    bool operator==(const Compound_Selector& rhs) const;
    bool operator<(const Simple_Selector& rhs) const;
    bool operator==(const Simple_Selector& rhs) const;

    ComplexSelectorSet& sources() { return sources_; }
    void clearSources() { sources_.clear(); }
    void mergeSources(ComplexSelectorSet& sources);

    Compound_Selector_Ptr minus(Compound_Selector_Ptr rhs);
    void cloneChildren() override;
    ATTACH_AST_OPERATIONS(Compound_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ////////////////////////////////////////////////////////////////////////////
  // General selectors -- i.e., simple sequences combined with one of the four
  // CSS selector combinators (">", "+", "~", and whitespace). Essentially a
  // linked list.
  ////////////////////////////////////////////////////////////////////////////
  class Complex_Selector final : public Selector {
  public:
    enum Combinator { ANCESTOR_OF, PARENT_OF, PRECEDES, ADJACENT_TO, REFERENCE };
  private:
    HASH_CONSTREF(Combinator, combinator)
    HASH_PROPERTY(Compound_Selector_Obj, head)
    HASH_PROPERTY(Complex_Selector_Obj, tail)
    HASH_PROPERTY(String_Obj, reference);
  public:
    bool contains_placeholder() {
      if (head() && head()->contains_placeholder()) return true;
      if (tail() && tail()->contains_placeholder()) return true;
      return false;
    };
    Complex_Selector(ParserState pstate,
                     Combinator c = ANCESTOR_OF,
                     Compound_Selector_Obj h = {},
                     Complex_Selector_Obj t = {},
                     String_Obj r = {});

    bool empty() const;

    bool has_parent_ref() const override;
    bool has_real_parent_ref() const override;
    Complex_Selector_Obj skip_empty_reference();

    // can still have a tail
    bool is_empty_ancestor() const;

    Selector_List_Ptr tails(Selector_List_Ptr tails);

    // front returns the first real tail
    // skips over parent and empty ones
    Complex_Selector_Ptr_Const first() const;
    Complex_Selector_Ptr mutable_first();

    // last returns the last real tail
    Complex_Selector_Ptr_Const last() const;
    Complex_Selector_Ptr mutable_last();

    size_t length() const;
    Selector_List_Ptr resolve_parent_refs(SelectorStack& pstack, Backtraces& traces, bool implicit_parent = true);
    bool is_superselector_of(Compound_Selector_Ptr_Const sub, std::string wrapping = "") const;
    bool is_superselector_of(Complex_Selector_Ptr_Const sub, std::string wrapping = "") const;
    bool is_superselector_of(Selector_List_Ptr_Const sub, std::string wrapping = "") const;
    Selector_List_Ptr unify_with(Complex_Selector_Ptr rhs);
    Combinator clear_innermost();
    void append(Complex_Selector_Obj, Backtraces& traces);
    void set_innermost(Complex_Selector_Obj, Combinator);

    size_t hash() const override;
    virtual unsigned long specificity() const override;
    virtual void set_media_block(Media_Block_Ptr mb) override;
    virtual bool has_placeholder();
    int unification_order() const override
    {
      throw std::runtime_error("unification_order for Complex_Selector is undefined");
    }
    bool find ( bool (*f)(AST_Node_Obj) ) override;

    bool operator<(const Selector& rhs) const override;
    bool operator==(const Selector& rhs) const override;
    bool operator<(const Selector_List& rhs) const;
    bool operator==(const Selector_List& rhs) const;
    bool operator<(const Complex_Selector& rhs) const;
    bool operator==(const Complex_Selector& rhs) const;
    bool operator<(const Compound_Selector& rhs) const;
    bool operator==(const Compound_Selector& rhs) const;
    bool operator<(const Simple_Selector& rhs) const;
    bool operator==(const Simple_Selector& rhs) const;

    const ComplexSelectorSet sources();
    void addSources(ComplexSelectorSet& sources);
    void clearSources();

    void cloneChildren() override;
    ATTACH_AST_OPERATIONS(Complex_Selector)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  ///////////////////////////////////
  // Comma-separated selector groups.
  ///////////////////////////////////
  class Selector_List final : public Selector, public Vectorized<Complex_Selector_Obj> {
    ADD_PROPERTY(Selector_Schema_Obj, schema)
    ADD_CONSTREF(std::vector<std::string>, wspace)
  protected:
    void adjust_after_pushing(Complex_Selector_Obj c) override;
  public:
    Selector_List(ParserState pstate, size_t s = 0);
    std::string type() const override { return "list"; }
    // remove parent selector references
    // basically unwraps parsed selectors
    bool has_parent_ref() const override;
    bool has_real_parent_ref() const override;
    void remove_parent_selectors();
    Selector_List_Ptr resolve_parent_refs(SelectorStack& pstack, Backtraces& traces, bool implicit_parent = true);
    bool is_superselector_of(Compound_Selector_Ptr_Const sub, std::string wrapping = "") const;
    bool is_superselector_of(Complex_Selector_Ptr_Const sub, std::string wrapping = "") const;
    bool is_superselector_of(Selector_List_Ptr_Const sub, std::string wrapping = "") const;
    Selector_List_Ptr unify_with(Selector_List_Ptr);
    void populate_extends(Selector_List_Obj, Subset_Map&);
    Selector_List_Obj eval(Eval& eval);

    size_t hash() const override;
    virtual unsigned long specificity() const override;
    virtual void set_media_block(Media_Block_Ptr mb) override;
    virtual bool has_placeholder();
    int unification_order() const override
    {
      throw std::runtime_error("unification_order for Selector_List is undefined");
    }
    bool find ( bool (*f)(AST_Node_Obj) ) override;
    bool operator<(const Selector& rhs) const override;
    bool operator==(const Selector& rhs) const override;
    bool operator<(const Selector_List& rhs) const;
    bool operator==(const Selector_List& rhs) const;
    bool operator<(const Complex_Selector& rhs) const;
    bool operator==(const Complex_Selector& rhs) const;
    bool operator<(const Compound_Selector& rhs) const;
    bool operator==(const Compound_Selector& rhs) const;
    bool operator<(const Simple_Selector& rhs) const;
    bool operator==(const Simple_Selector& rhs) const;
    // Selector Lists can be compared to comma lists
    bool operator<(const Expression& rhs) const override;
    bool operator==(const Expression& rhs) const override;
    void cloneChildren() override;
    ATTACH_AST_OPERATIONS(Selector_List)
    ATTACH_CRTP_PERFORM_METHODS()
  };

  // compare function for sorting and probably other other uses
  struct cmp_complex_selector { inline bool operator() (const Complex_Selector_Obj l, const Complex_Selector_Obj r) { return (*l < *r); } };
  struct cmp_compound_selector { inline bool operator() (const Compound_Selector_Obj l, const Compound_Selector_Obj r) { return (*l < *r); } };
  struct cmp_simple_selector { inline bool operator() (const Simple_Selector_Obj l, const Simple_Selector_Obj r) { return (*l < *r); } };

}

#endif

#ifdef __clang__

// #pragma clang diagnostic pop
// #pragma clang diagnostic push

#endif

#endif
#ifndef SASS_NODE_H
#define SASS_NODE_H

#include <memory>



namespace Sass {




  class Context;

  /*
   There are a lot of stumbling blocks when trying to port the ruby extend code to C++. The biggest is the choice of
   data type. The ruby code will pretty seamlessly switch types between an Array<SimpleSequence or Op> (libsass'
   equivalent is the Complex_Selector) to a Sequence, which contains more metadata about the sequence than just the
   selector info. They also have the ability to have arbitrary nestings of arrays like [1, [2]], which is hard to
   implement using Array equivalents in C++ (like the deque or vector). They also have the ability to include nil
   in the arrays, like [1, nil, 3], which has potential semantic differences than an empty array [1, [], 3]. To be
   able to represent all of these as unique cases, we need to create a tree of variant objects. The tree nature allows
   the inconsistent nesting levels. The variant nature (while making some of the C++ code uglier) allows the code to
   more closely match the ruby code, which is a huge benefit when attempting to implement an complex algorithm like
   the Extend operator.

   Note that the current libsass data model also pairs the combinator with the Complex_Selector that follows it, but
   ruby sass has no such restriction, so we attempt to create a data structure that can handle them split apart.
   */

  class Node;
  typedef std::deque<Node> NodeDeque;
  typedef std::shared_ptr<NodeDeque> NodeDequePtr;

  class Node {
  public:
    enum TYPE {
      SELECTOR,
      COMBINATOR,
      COLLECTION,
      NIL
    };

    TYPE type() const { return mType; }
    bool isCombinator() const { return mType == COMBINATOR; }
    bool isSelector() const { return mType == SELECTOR; }
    bool isCollection() const { return mType == COLLECTION; }
    bool isNil() const { return mType == NIL; }
    bool got_line_feed;

    Complex_Selector::Combinator combinator() const { return mCombinator; }

    Complex_Selector_Obj selector() { return mpSelector; }
    Complex_Selector_Obj selector() const { return mpSelector; }

    NodeDequePtr collection() { return mpCollection; }
    const NodeDequePtr collection() const { return mpCollection; }

    static Node createCombinator(const Complex_Selector::Combinator& combinator);

    // This method will klone the selector, stripping off the tail and combinator
    static Node createSelector(const Complex_Selector& pSelector);

    static Node createCollection();
    static Node createCollection(const NodeDeque& values);

    static Node createNil();
    static Node naiveTrim(Node& seqses);

    Node klone() const;

    bool operator==(const Node& rhs) const;
    inline bool operator!=(const Node& rhs) const { return !(*this == rhs); }


    /*
    COLLECTION FUNCTIONS

    Most types don't need any helper methods (nil and combinator due to their simplicity and
    selector due to the fact that we leverage the non-node selector code on the Complex_Selector
    whereever possible). The following methods are intended to be called on Node objects whose
    type is COLLECTION only.
    */

    // rhs and this must be node collections. Shallow copy the nodes from rhs to the end of this.
    // This function DOES NOT remove the nodes from rhs.
    void plus(Node& rhs);

    // potentialChild must be a node collection of selectors/combinators. this must be a collection
    // of collections of nodes/combinators. This method checks if potentialChild is a child of this
    // Node.
    bool contains(const Node& potentialChild) const;

  private:
    // Private constructor; Use the static methods (like createCombinator and createSelector)
    // to instantiate this object. This is more expressive, and it allows us to break apart each
    // case into separate functions.
    Node(const TYPE& type, Complex_Selector::Combinator combinator, Complex_Selector_Ptr pSelector, NodeDequePtr& pCollection);

    TYPE mType;

    // TODO: can we union these to save on memory?
    Complex_Selector::Combinator mCombinator;
    Complex_Selector_Obj mpSelector;
    NodeDequePtr mpCollection;
  };

#ifdef DEBUG
  std::ostream& operator<<(std::ostream& os, const Node& node);
#endif
  Node complexSelectorToNode(Complex_Selector_Ptr pToConvert);
  Complex_Selector_Ptr nodeToComplexSelector(const Node& toConvert);

}

#endif
#ifndef SASS_EVAL_H
#define SASS_EVAL_H

#ifndef SASS_LISTIZE_H
#define SASS_LISTIZE_H

// sass.hpp must go before all system headers to get the
// __EXTENSIONS__ fix on Solaris.



namespace Sass {

  struct Backtrace;

  class Listize : public Operation_CRTP<Expression_Ptr, Listize> {

  public:
    Listize();
    ~Listize() { }

    Expression_Ptr operator()(Selector_List_Ptr);
    Expression_Ptr operator()(Complex_Selector_Ptr);
    Expression_Ptr operator()(Compound_Selector_Ptr);

    // generic fallback
    template <typename U>
    Expression_Ptr fallback(U x)
    { return Cast<Expression>(x); }
  };

}

#endif

namespace Sass {

  class Expand;
  class Context;

  class Eval : public Operation_CRTP<Expression_Ptr, Eval> {

   public:
    Expand& exp;
    Context& ctx;
    Backtraces& traces;
    Eval(Expand& exp);
    ~Eval();

    bool force;
    bool is_in_comment;
    bool is_in_selector_schema;

    Boolean_Obj bool_true;
    Boolean_Obj bool_false;

    Env* environment();
    EnvStack& env_stack();
    const std::string cwd();
    Selector_List_Obj selector();
    CalleeStack& callee_stack();
    SelectorStack& selector_stack();
    bool& old_at_root_without_rule();
    struct Sass_Inspect_Options& options();
    struct Sass_Inspect_Options options2();
    struct Sass_Compiler* compiler();

    // for evaluating function bodies
    Expression_Ptr operator()(Block_Ptr);
    Expression_Ptr operator()(Assignment_Ptr);
    Expression_Ptr operator()(If_Ptr);
    Expression_Ptr operator()(For_Ptr);
    Expression_Ptr operator()(Each_Ptr);
    Expression_Ptr operator()(While_Ptr);
    Expression_Ptr operator()(Return_Ptr);
    Expression_Ptr operator()(Warning_Ptr);
    Expression_Ptr operator()(Error_Ptr);
    Expression_Ptr operator()(Debug_Ptr);

    Expression_Ptr operator()(List_Ptr);
    Expression_Ptr operator()(Map_Ptr);
    Expression_Ptr operator()(Binary_Expression_Ptr);
    Expression_Ptr operator()(Unary_Expression_Ptr);
    Expression_Ptr operator()(Function_Call_Ptr);
    Expression_Ptr operator()(Variable_Ptr);
    Expression_Ptr operator()(Number_Ptr);
    Expression_Ptr operator()(Color_RGBA_Ptr);
    Expression_Ptr operator()(Color_HSLA_Ptr);
    Expression_Ptr operator()(Boolean_Ptr);
    Expression_Ptr operator()(String_Schema_Ptr);
    Expression_Ptr operator()(String_Quoted_Ptr);
    Expression_Ptr operator()(String_Constant_Ptr);
    // Expression_Ptr operator()(Selector_List_Ptr);
    Media_Query_Ptr operator()(Media_Query_Ptr);
    Expression_Ptr operator()(Media_Query_Expression_Ptr);
    Expression_Ptr operator()(At_Root_Query_Ptr);
    Expression_Ptr operator()(Supports_Operator_Ptr);
    Expression_Ptr operator()(Supports_Negation_Ptr);
    Expression_Ptr operator()(Supports_Declaration_Ptr);
    Expression_Ptr operator()(Supports_Interpolation_Ptr);
    Expression_Ptr operator()(Null_Ptr);
    Expression_Ptr operator()(Argument_Ptr);
    Expression_Ptr operator()(Arguments_Ptr);
    Expression_Ptr operator()(Comment_Ptr);

    // these will return selectors
    Selector_List_Ptr operator()(Selector_List_Ptr);
    Selector_List_Ptr operator()(Complex_Selector_Ptr);
    Compound_Selector_Ptr operator()(Compound_Selector_Ptr);
    Simple_Selector_Ptr operator()(Simple_Selector_Ptr s);
    Wrapped_Selector_Ptr operator()(Wrapped_Selector_Ptr s);

    // they don't have any specific implementation (yet)
    Id_Selector_Ptr operator()(Id_Selector_Ptr s) { return s; };
    Class_Selector_Ptr operator()(Class_Selector_Ptr s) { return s; };
    Pseudo_Selector_Ptr operator()(Pseudo_Selector_Ptr s) { return s; };
    Type_Selector_Ptr operator()(Type_Selector_Ptr s) { return s; };
    Attribute_Selector_Ptr operator()(Attribute_Selector_Ptr s) { return s; };
    Placeholder_Selector_Ptr operator()(Placeholder_Selector_Ptr s) { return s; };

    // actual evaluated selectors
    Selector_List_Ptr operator()(Selector_Schema_Ptr);
    Expression_Ptr operator()(Parent_Selector_Ptr);
    Expression_Ptr operator()(Parent_Reference_Ptr);

    // generic fallback
    template <typename U>
    Expression_Ptr fallback(U x)
    { return Cast<Expression>(x); }

  private:
    void interpolation(Context& ctx, std::string& res, Expression_Obj ex, bool into_quotes, bool was_itpl = false);

  };

}

#endif

namespace Sass {

  Node subweave(Node& one, Node& two);

  class Extend : public Operation_CRTP<void, Extend> {

    Subset_Map& subset_map;
    Eval* eval;

  private:

    std::unordered_map<
      Selector_List_Obj, // key
      Selector_List_Obj, // value
      HashNodes, // hasher
      CompareNodes // compare
    > memoizeList;

    std::unordered_map<
      Complex_Selector_Obj, // key
      Node, // value
      HashNodes, // hasher
      CompareNodes // compare
    > memoizeComplex;

    /* this turned out to be too much overhead
       re-evaluate once we store an ast selector
    std::unordered_map<
      Compound_Selector_Obj, // key
      Node, // value
      HashNodes, // hasher
      CompareNodes // compare
    > memoizeCompound;
    */

    void extendObjectWithSelectorAndBlock(Ruleset_Ptr pObject);
    Node extendComplexSelector(Complex_Selector_Ptr sel, CompoundSelectorSet& seen, bool isReplace, bool isOriginal);
    Node extendCompoundSelector(Compound_Selector_Ptr sel, CompoundSelectorSet& seen, bool isReplace);
    bool complexSelectorHasExtension(Complex_Selector_Ptr selector, CompoundSelectorSet& seen);
    Node trim(Node& seqses, bool isReplace);
    Node weave(Node& path);

  public:
    void setEval(Eval& eval);
    Selector_List_Ptr extendSelectorList(Selector_List_Obj pSelectorList, bool isReplace, bool& extendedSomething, CompoundSelectorSet& seen);
    Selector_List_Ptr extendSelectorList(Selector_List_Obj pSelectorList, bool isReplace = false) {
      bool extendedSomething = false;
      CompoundSelectorSet seen;
      return extendSelectorList(pSelectorList, isReplace, extendedSomething, seen);
    }
    Selector_List_Ptr extendSelectorList(Selector_List_Obj pSelectorList, CompoundSelectorSet& seen) {
      bool isReplace = false;
      bool extendedSomething = false;
      return extendSelectorList(pSelectorList, isReplace, extendedSomething, seen);
    }
    Extend(Subset_Map&);
    ~Extend() { }

    void operator()(Block_Ptr);
    void operator()(Ruleset_Ptr);
    void operator()(Supports_Block_Ptr);
    void operator()(Media_Block_Ptr);
    void operator()(Directive_Ptr);

    // ignore missed types
    template <typename U>
    void fallback(U x) {}

  };

}

#endif
#ifndef SASS_OPERATORS_H
#define SASS_OPERATORS_H

#ifndef SASS_VALUES_H
#define SASS_VALUES_H


namespace Sass {

  union Sass_Value* ast_node_to_sass_value (const Expression_Ptr val);
  Value_Ptr sass_value_to_ast_node (const union Sass_Value* val);

}
#endif

namespace Sass {

  namespace Operators {

    // equality operator using AST Node operator==
    bool eq(Expression_Obj, Expression_Obj);
    bool neq(Expression_Obj, Expression_Obj);
    // specific operators based on cmp and eq
    bool lt(Expression_Obj, Expression_Obj);
    bool gt(Expression_Obj, Expression_Obj);
    bool lte(Expression_Obj, Expression_Obj);
    bool gte(Expression_Obj, Expression_Obj);
    // arithmetic for all the combinations that matter
    Value_Ptr op_strings(Sass::Operand, Value&, Value&, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed = false);
    Value_Ptr op_colors(enum Sass_OP, const Color_RGBA&, const Color_RGBA&, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed = false);
    Value_Ptr op_numbers(enum Sass_OP, const Number&, const Number&, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed = false);
    Value_Ptr op_number_color(enum Sass_OP, const Number&, const Color_RGBA&, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed = false);
    Value_Ptr op_color_number(enum Sass_OP, const Color_RGBA&, const Number&, struct Sass_Inspect_Options opt, const ParserState& pstate, bool delayed = false);

  };

}

#endif
#ifndef SASS_FN_COLORS_H
#define SASS_FN_COLORS_H


namespace Sass {

  namespace Functions {

    // macros for common ranges (u mean unsigned or upper, r for full range)
    #define DARG_U_FACT(argname) get_arg_r(argname, env, sig, pstate, traces, - 0.0, 1.0) // double
    #define DARG_R_FACT(argname) get_arg_r(argname, env, sig, pstate, traces, - 1.0, 1.0) // double
    #define DARG_U_BYTE(argname) get_arg_r(argname, env, sig, pstate, traces, - 0.0, 255.0) // double
    #define DARG_R_BYTE(argname) get_arg_r(argname, env, sig, pstate, traces, - 255.0, 255.0) // double
    #define DARG_U_PRCT(argname) get_arg_r(argname, env, sig, pstate, traces, - 0.0, 100.0) // double
    #define DARG_R_PRCT(argname) get_arg_r(argname, env, sig, pstate, traces, - 100.0, 100.0) // double

    // macros for color related inputs (rbg and alpha/opacity values)
    #define COLOR_NUM(argname) color_num(argname, env, sig, pstate, traces) // double
    #define ALPHA_NUM(argname) alpha_num(argname, env, sig, pstate, traces) // double

    extern Signature rgb_sig;
    extern Signature rgba_4_sig;
    extern Signature rgba_2_sig;
    extern Signature red_sig;
    extern Signature green_sig;
    extern Signature blue_sig;
    extern Signature mix_sig;
    extern Signature hsl_sig;
    extern Signature hsla_sig;
    extern Signature hue_sig;
    extern Signature saturation_sig;
    extern Signature lightness_sig;
    extern Signature adjust_hue_sig;
    extern Signature lighten_sig;
    extern Signature darken_sig;
    extern Signature saturate_sig;
    extern Signature desaturate_sig;
    extern Signature grayscale_sig;
    extern Signature complement_sig;
    extern Signature invert_sig;
    extern Signature alpha_sig;
    extern Signature opacity_sig;
    extern Signature opacify_sig;
    extern Signature fade_in_sig;
    extern Signature transparentize_sig;
    extern Signature fade_out_sig;
    extern Signature adjust_color_sig;
    extern Signature scale_color_sig;
    extern Signature change_color_sig;
    extern Signature ie_hex_str_sig;

    BUILT_IN(rgb);
    BUILT_IN(rgba_4);
    BUILT_IN(rgba_2);
    BUILT_IN(red);
    BUILT_IN(green);
    BUILT_IN(blue);
    BUILT_IN(mix);
    BUILT_IN(hsl);
    BUILT_IN(hsla);
    BUILT_IN(hue);
    BUILT_IN(saturation);
    BUILT_IN(lightness);
    BUILT_IN(adjust_hue);
    BUILT_IN(lighten);
    BUILT_IN(darken);
    BUILT_IN(saturate);
    BUILT_IN(desaturate);
    BUILT_IN(grayscale);
    BUILT_IN(complement);
    BUILT_IN(invert);
    BUILT_IN(alpha);
    BUILT_IN(opacify);
    BUILT_IN(transparentize);
    BUILT_IN(adjust_color);
    BUILT_IN(scale_color);
    BUILT_IN(change_color);
    BUILT_IN(ie_hex_str);

  }

}

#endif
#ifndef SASS_DEBUG_H
#define SASS_DEBUG_H

#include <stdint.h>

#ifndef UINT32_MAX
  #define UINT32_MAX 0xffffffffU
#endif

enum dbg_lvl_t : uint32_t {
  NONE = 0,
  TRIM = 1,
  CHUNKS = 2,
  SUBWEAVE = 4,
  WEAVE = 8,
  EXTEND_COMPOUND = 16,
  EXTEND_COMPLEX = 32,
  LCS = 64,
  EXTEND_OBJECT = 128,
  ALL = UINT32_MAX
};

#ifdef DEBUG

#ifndef DEBUG_LVL
const uint32_t debug_lvl = UINT32_MAX;
#else
const uint32_t debug_lvl = (DEBUG_LVL);
#endif // DEBUG_LVL

#define DEBUG_PRINT(lvl, x) if((lvl) & debug_lvl) { std::cerr << x; }
#define DEBUG_PRINTLN(lvl, x) if((lvl) & debug_lvl) { std::cerr << x << std::endl; }
#define DEBUG_EXEC(lvl, x) if((lvl) & debug_lvl) { x; }

#else // DEBUG

#define DEBUG_PRINT(lvl, x)
#define DEBUG_PRINTLN(lvl, x)
#define DEBUG_EXEC(lvl, x)

#endif // DEBUG

#endif // SASS_DEBUG

#ifndef SASS_COLOR_MAPS_H
#define SASS_COLOR_MAPS_H


namespace Sass {

  struct map_cmp_str
  {
    bool operator()(char const *a, char const *b) const
    {
      return std::strcmp(a, b) < 0;
    }
  };

  namespace ColorNames
  {
    extern const char aliceblue[];
    extern const char antiquewhite[];
    extern const char cyan[];
    extern const char aqua[];
    extern const char aquamarine[];
    extern const char azure[];
    extern const char beige[];
    extern const char bisque[];
    extern const char black[];
    extern const char blanchedalmond[];
    extern const char blue[];
    extern const char blueviolet[];
    extern const char brown[];
    extern const char burlywood[];
    extern const char cadetblue[];
    extern const char chartreuse[];
    extern const char chocolate[];
    extern const char coral[];
    extern const char cornflowerblue[];
    extern const char cornsilk[];
    extern const char crimson[];
    extern const char darkblue[];
    extern const char darkcyan[];
    extern const char darkgoldenrod[];
    extern const char darkgray[];
    extern const char darkgrey[];
    extern const char darkgreen[];
    extern const char darkkhaki[];
    extern const char darkmagenta[];
    extern const char darkolivegreen[];
    extern const char darkorange[];
    extern const char darkorchid[];
    extern const char darkred[];
    extern const char darksalmon[];
    extern const char darkseagreen[];
    extern const char darkslateblue[];
    extern const char darkslategray[];
    extern const char darkslategrey[];
    extern const char darkturquoise[];
    extern const char darkviolet[];
    extern const char deeppink[];
    extern const char deepskyblue[];
    extern const char dimgray[];
    extern const char dimgrey[];
    extern const char dodgerblue[];
    extern const char firebrick[];
    extern const char floralwhite[];
    extern const char forestgreen[];
    extern const char magenta[];
    extern const char fuchsia[];
    extern const char gainsboro[];
    extern const char ghostwhite[];
    extern const char gold[];
    extern const char goldenrod[];
    extern const char gray[];
    extern const char grey[];
    extern const char green[];
    extern const char greenyellow[];
    extern const char honeydew[];
    extern const char hotpink[];
    extern const char indianred[];
    extern const char indigo[];
    extern const char ivory[];
    extern const char khaki[];
    extern const char lavender[];
    extern const char lavenderblush[];
    extern const char lawngreen[];
    extern const char lemonchiffon[];
    extern const char lightblue[];
    extern const char lightcoral[];
    extern const char lightcyan[];
    extern const char lightgoldenrodyellow[];
    extern const char lightgray[];
    extern const char lightgrey[];
    extern const char lightgreen[];
    extern const char lightpink[];
    extern const char lightsalmon[];
    extern const char lightseagreen[];
    extern const char lightskyblue[];
    extern const char lightslategray[];
    extern const char lightslategrey[];
    extern const char lightsteelblue[];
    extern const char lightyellow[];
    extern const char lime[];
    extern const char limegreen[];
    extern const char linen[];
    extern const char maroon[];
    extern const char mediumaquamarine[];
    extern const char mediumblue[];
    extern const char mediumorchid[];
    extern const char mediumpurple[];
    extern const char mediumseagreen[];
    extern const char mediumslateblue[];
    extern const char mediumspringgreen[];
    extern const char mediumturquoise[];
    extern const char mediumvioletred[];
    extern const char midnightblue[];
    extern const char mintcream[];
    extern const char mistyrose[];
    extern const char moccasin[];
    extern const char navajowhite[];
    extern const char navy[];
    extern const char oldlace[];
    extern const char olive[];
    extern const char olivedrab[];
    extern const char orange[];
    extern const char orangered[];
    extern const char orchid[];
    extern const char palegoldenrod[];
    extern const char palegreen[];
    extern const char paleturquoise[];
    extern const char palevioletred[];
    extern const char papayawhip[];
    extern const char peachpuff[];
    extern const char peru[];
    extern const char pink[];
    extern const char plum[];
    extern const char powderblue[];
    extern const char purple[];
    extern const char red[];
    extern const char rosybrown[];
    extern const char royalblue[];
    extern const char saddlebrown[];
    extern const char salmon[];
    extern const char sandybrown[];
    extern const char seagreen[];
    extern const char seashell[];
    extern const char sienna[];
    extern const char silver[];
    extern const char skyblue[];
    extern const char slateblue[];
    extern const char slategray[];
    extern const char slategrey[];
    extern const char snow[];
    extern const char springgreen[];
    extern const char steelblue[];
    extern const char tan[];
    extern const char teal[];
    extern const char thistle[];
    extern const char tomato[];
    extern const char turquoise[];
    extern const char violet[];
    extern const char wheat[];
    extern const char white[];
    extern const char whitesmoke[];
    extern const char yellow[];
    extern const char yellowgreen[];
    extern const char rebeccapurple[];
    extern const char transparent[];
  }

  namespace Colors {
    extern const Color_RGBA aliceblue;
    extern const Color_RGBA antiquewhite;
    extern const Color_RGBA cyan;
    extern const Color_RGBA aqua;
    extern const Color_RGBA aquamarine;
    extern const Color_RGBA azure;
    extern const Color_RGBA beige;
    extern const Color_RGBA bisque;
    extern const Color_RGBA black;
    extern const Color_RGBA blanchedalmond;
    extern const Color_RGBA blue;
    extern const Color_RGBA blueviolet;
    extern const Color_RGBA brown;
    extern const Color_RGBA burlywood;
    extern const Color_RGBA cadetblue;
    extern const Color_RGBA chartreuse;
    extern const Color_RGBA chocolate;
    extern const Color_RGBA coral;
    extern const Color_RGBA cornflowerblue;
    extern const Color_RGBA cornsilk;
    extern const Color_RGBA crimson;
    extern const Color_RGBA darkblue;
    extern const Color_RGBA darkcyan;
    extern const Color_RGBA darkgoldenrod;
    extern const Color_RGBA darkgray;
    extern const Color_RGBA darkgrey;
    extern const Color_RGBA darkgreen;
    extern const Color_RGBA darkkhaki;
    extern const Color_RGBA darkmagenta;
    extern const Color_RGBA darkolivegreen;
    extern const Color_RGBA darkorange;
    extern const Color_RGBA darkorchid;
    extern const Color_RGBA darkred;
    extern const Color_RGBA darksalmon;
    extern const Color_RGBA darkseagreen;
    extern const Color_RGBA darkslateblue;
    extern const Color_RGBA darkslategray;
    extern const Color_RGBA darkslategrey;
    extern const Color_RGBA darkturquoise;
    extern const Color_RGBA darkviolet;
    extern const Color_RGBA deeppink;
    extern const Color_RGBA deepskyblue;
    extern const Color_RGBA dimgray;
    extern const Color_RGBA dimgrey;
    extern const Color_RGBA dodgerblue;
    extern const Color_RGBA firebrick;
    extern const Color_RGBA floralwhite;
    extern const Color_RGBA forestgreen;
    extern const Color_RGBA magenta;
    extern const Color_RGBA fuchsia;
    extern const Color_RGBA gainsboro;
    extern const Color_RGBA ghostwhite;
    extern const Color_RGBA gold;
    extern const Color_RGBA goldenrod;
    extern const Color_RGBA gray;
    extern const Color_RGBA grey;
    extern const Color_RGBA green;
    extern const Color_RGBA greenyellow;
    extern const Color_RGBA honeydew;
    extern const Color_RGBA hotpink;
    extern const Color_RGBA indianred;
    extern const Color_RGBA indigo;
    extern const Color_RGBA ivory;
    extern const Color_RGBA khaki;
    extern const Color_RGBA lavender;
    extern const Color_RGBA lavenderblush;
    extern const Color_RGBA lawngreen;
    extern const Color_RGBA lemonchiffon;
    extern const Color_RGBA lightblue;
    extern const Color_RGBA lightcoral;
    extern const Color_RGBA lightcyan;
    extern const Color_RGBA lightgoldenrodyellow;
    extern const Color_RGBA lightgray;
    extern const Color_RGBA lightgrey;
    extern const Color_RGBA lightgreen;
    extern const Color_RGBA lightpink;
    extern const Color_RGBA lightsalmon;
    extern const Color_RGBA lightseagreen;
    extern const Color_RGBA lightskyblue;
    extern const Color_RGBA lightslategray;
    extern const Color_RGBA lightslategrey;
    extern const Color_RGBA lightsteelblue;
    extern const Color_RGBA lightyellow;
    extern const Color_RGBA lime;
    extern const Color_RGBA limegreen;
    extern const Color_RGBA linen;
    extern const Color_RGBA maroon;
    extern const Color_RGBA mediumaquamarine;
    extern const Color_RGBA mediumblue;
    extern const Color_RGBA mediumorchid;
    extern const Color_RGBA mediumpurple;
    extern const Color_RGBA mediumseagreen;
    extern const Color_RGBA mediumslateblue;
    extern const Color_RGBA mediumspringgreen;
    extern const Color_RGBA mediumturquoise;
    extern const Color_RGBA mediumvioletred;
    extern const Color_RGBA midnightblue;
    extern const Color_RGBA mintcream;
    extern const Color_RGBA mistyrose;
    extern const Color_RGBA moccasin;
    extern const Color_RGBA navajowhite;
    extern const Color_RGBA navy;
    extern const Color_RGBA oldlace;
    extern const Color_RGBA olive;
    extern const Color_RGBA olivedrab;
    extern const Color_RGBA orange;
    extern const Color_RGBA orangered;
    extern const Color_RGBA orchid;
    extern const Color_RGBA palegoldenrod;
    extern const Color_RGBA palegreen;
    extern const Color_RGBA paleturquoise;
    extern const Color_RGBA palevioletred;
    extern const Color_RGBA papayawhip;
    extern const Color_RGBA peachpuff;
    extern const Color_RGBA peru;
    extern const Color_RGBA pink;
    extern const Color_RGBA plum;
    extern const Color_RGBA powderblue;
    extern const Color_RGBA purple;
    extern const Color_RGBA red;
    extern const Color_RGBA rosybrown;
    extern const Color_RGBA royalblue;
    extern const Color_RGBA saddlebrown;
    extern const Color_RGBA salmon;
    extern const Color_RGBA sandybrown;
    extern const Color_RGBA seagreen;
    extern const Color_RGBA seashell;
    extern const Color_RGBA sienna;
    extern const Color_RGBA silver;
    extern const Color_RGBA skyblue;
    extern const Color_RGBA slateblue;
    extern const Color_RGBA slategray;
    extern const Color_RGBA slategrey;
    extern const Color_RGBA snow;
    extern const Color_RGBA springgreen;
    extern const Color_RGBA steelblue;
    extern const Color_RGBA tan;
    extern const Color_RGBA teal;
    extern const Color_RGBA thistle;
    extern const Color_RGBA tomato;
    extern const Color_RGBA turquoise;
    extern const Color_RGBA violet;
    extern const Color_RGBA wheat;
    extern const Color_RGBA white;
    extern const Color_RGBA whitesmoke;
    extern const Color_RGBA yellow;
    extern const Color_RGBA yellowgreen;
    extern const Color_RGBA rebeccapurple;
    extern const Color_RGBA transparent;
  }

  Color_RGBA_Ptr_Const name_to_color(const char*);
  Color_RGBA_Ptr_Const name_to_color(const std::string&);
  const char* color_to_name(const int);
  const char* color_to_name(const Color_RGBA&);
  const char* color_to_name(const double);

}

#endif
#ifndef SASS_SASS_FUNCTIONS_H
#define SASS_SASS_FUNCTIONS_H


// Struct to hold custom function callback
struct Sass_Function {
  char*            signature;
  Sass_Function_Fn function;
  void*            cookie;
};

// External import entry
struct Sass_Import {
  char* imp_path; // path as found in the import statement
  char *abs_path; // path after importer has resolved it
  char* source;
  char* srcmap;
  // error handling
  char* error;
  size_t line;
  size_t column;
};

// External environments
struct Sass_Env {
  // links to parent frames
  Sass::Env* frame;
};

// External call entry
struct Sass_Callee {
  const char* name;
  const char* path;
  size_t line;
  size_t column;
  enum Sass_Callee_Type type;
  struct Sass_Env env;
};

// Struct to hold importer callback
struct Sass_Importer {
  Sass_Importer_Fn importer;
  double           priority;
  void*            cookie;
};

#endif
#ifndef SASS_TO_VALUE_H
#define SASS_TO_VALUE_H


namespace Sass {

  class To_Value : public Operation_CRTP<Value_Ptr, To_Value> {

  private:

    Context& ctx;

  public:

    To_Value(Context& ctx)
    : ctx(ctx)
    { }
    ~To_Value() { }
    using Operation<Value_Ptr>::operator();

    Value_Ptr operator()(Argument_Ptr);
    Value_Ptr operator()(Boolean_Ptr);
    Value_Ptr operator()(Number_Ptr);
    Value_Ptr operator()(Color_RGBA_Ptr);
    Value_Ptr operator()(Color_HSLA_Ptr);
    Value_Ptr operator()(String_Constant_Ptr);
    Value_Ptr operator()(String_Quoted_Ptr);
    Value_Ptr operator()(Custom_Warning_Ptr);
    Value_Ptr operator()(Custom_Error_Ptr);
    Value_Ptr operator()(List_Ptr);
    Value_Ptr operator()(Map_Ptr);
    Value_Ptr operator()(Null_Ptr);
    Value_Ptr operator()(Function_Ptr);

    // convert to string via `To_String`
    Value_Ptr operator()(Selector_List_Ptr);
    Value_Ptr operator()(Binary_Expression_Ptr);

  };

}

#endif
#ifndef SASS_EXPAND_H
#define SASS_EXPAND_H



namespace Sass {

  class Listize;
  class Context;
  class Eval;
  struct Backtrace;

  class Expand : public Operation_CRTP<Statement_Ptr, Expand> {
  public:

    Env* environment();
    Selector_List_Obj selector();

    Context&          ctx;
    Backtraces&       traces;
    Eval              eval;
    size_t            recursions;
    bool              in_keyframes;
    bool              at_root_without_rule;
    bool              old_at_root_without_rule;

    // it's easier to work with vectors
    EnvStack      env_stack;
    BlockStack    block_stack;
    CallStack     call_stack;
    SelectorStack selector_stack;
    MediaStack    media_stack;

    Boolean_Obj bool_true;

  private:
    void expand_selector_list(Selector_Obj, Selector_List_Obj extender);

  public:
    Expand(Context&, Env*, SelectorStack* stack = NULL);
    ~Expand() { }

    Block_Ptr operator()(Block_Ptr);
    Statement_Ptr operator()(Ruleset_Ptr);
    Statement_Ptr operator()(Media_Block_Ptr);
    Statement_Ptr operator()(Supports_Block_Ptr);
    Statement_Ptr operator()(At_Root_Block_Ptr);
    Statement_Ptr operator()(Directive_Ptr);
    Statement_Ptr operator()(Declaration_Ptr);
    Statement_Ptr operator()(Assignment_Ptr);
    Statement_Ptr operator()(Import_Ptr);
    Statement_Ptr operator()(Import_Stub_Ptr);
    Statement_Ptr operator()(Warning_Ptr);
    Statement_Ptr operator()(Error_Ptr);
    Statement_Ptr operator()(Debug_Ptr);
    Statement_Ptr operator()(Comment_Ptr);
    Statement_Ptr operator()(If_Ptr);
    Statement_Ptr operator()(For_Ptr);
    Statement_Ptr operator()(Each_Ptr);
    Statement_Ptr operator()(While_Ptr);
    Statement_Ptr operator()(Return_Ptr);
    Statement_Ptr operator()(Extension_Ptr);
    Statement_Ptr operator()(Definition_Ptr);
    Statement_Ptr operator()(Mixin_Call_Ptr);
    Statement_Ptr operator()(Content_Ptr);

    void append_block(Block_Ptr);

  };

}

#endif
#ifndef SASS_PARSER_H
#define SASS_PARSER_H



#ifndef SASS_PRELEXER_H
#define SASS_PRELEXER_H

#ifndef SASS_LEXER_H
#define SASS_LEXER_H


namespace Sass {
  namespace Prelexer {

    //####################################
    // BASIC CHARACTER MATCHERS
    //####################################

    // Match standard control chars
    const char* kwd_at(const char* src);
    const char* kwd_dot(const char* src);
    const char* kwd_comma(const char* src);
    const char* kwd_colon(const char* src);
    const char* kwd_star(const char* src);
    const char* kwd_plus(const char* src);
    const char* kwd_minus(const char* src);
    const char* kwd_slash(const char* src);

    //####################################
    // BASIC CLASS MATCHERS
    //####################################

    // These are locale independant
    bool is_space(const char& src);
    bool is_alpha(const char& src);
    bool is_punct(const char& src);
    bool is_digit(const char& src);
    bool is_number(const char& src);
    bool is_alnum(const char& src);
    bool is_xdigit(const char& src);
    bool is_unicode(const char& src);
    bool is_nonascii(const char& src);
    bool is_character(const char& src);
    bool is_uri_character(const char& src);
    bool escapable_character(const char& src);

    // Match a single ctype predicate.
    const char* space(const char* src);
    const char* alpha(const char* src);
    const char* digit(const char* src);
    const char* xdigit(const char* src);
    const char* alnum(const char* src);
    const char* punct(const char* src);
    const char* hyphen(const char* src);
    const char* unicode(const char* src);
    const char* nonascii(const char* src);
    const char* character(const char* src);
    const char* uri_character(const char* src);
    const char* escapable_character(const char* src);

    // Match multiple ctype characters.
    const char* spaces(const char* src);
    const char* digits(const char* src);
    const char* hyphens(const char* src);

    // Whitespace handling.
    const char* no_spaces(const char* src);
    const char* optional_spaces(const char* src);

    // Match any single character (/./).
    const char* any_char(const char* src);

    // Assert word boundary (/\b/)
    // Is a zero-width positive lookaheads
    const char* word_boundary(const char* src);

    // Match a single linebreak (/(?:\n|\r\n?)/).
    const char* re_linebreak(const char* src);

    // Assert string boundaries (/\Z|\z|\A/)
    // There are zero-width positive lookaheads
    const char* end_of_line(const char* src);

    // Assert end_of_file boundary (/\z/)
    const char* end_of_file(const char* src);
    // const char* start_of_string(const char* src);

    // Type definition for prelexer functions
    typedef const char* (*prelexer)(const char*);

    //####################################
    // BASIC "REGEX" CONSTRUCTORS
    //####################################

    // Match a single character literal.
    // Regex equivalent: /(?:x)/
    template <char chr>
    const char* exactly(const char* src) {
      return *src == chr ? src + 1 : 0;
    }

    // Match the full string literal.
    // Regex equivalent: /(?:literal)/
    template <const char* str>
    const char* exactly(const char* src) {
      if (str == NULL) return 0;
      const char* pre = str;
      if (src == NULL) return 0;
      // there is a small chance that the search string
      // is longer than the rest of the string to look at
      while (*pre && *src == *pre) {
        ++src, ++pre;
      }
      // did the matcher finish?
      return *pre == 0 ? src : 0;
    }


    // Match a single character literal.
    // Regex equivalent: /(?:x)/i
    // only define lower case alpha chars
    template <char chr>
    const char* insensitive(const char* src) {
      return *src == chr || *src+32 == chr ? src + 1 : 0;
    }

    // Match the full string literal.
    // Regex equivalent: /(?:literal)/i
    // only define lower case alpha chars
    template <const char* str>
    const char* insensitive(const char* src) {
      if (str == NULL) return 0;
      const char* pre = str;
      if (src == NULL) return 0;
      // there is a small chance that the search string
      // is longer than the rest of the string to look at
      while (*pre && (*src == *pre || *src+32 == *pre)) {
        ++src, ++pre;
      }
      // did the matcher finish?
      return *pre == 0 ? src : 0;
    }

    // Match for members of char class.
    // Regex equivalent: /[axy]/
    template <const char* char_class>
    const char* class_char(const char* src) {
      const char* cc = char_class;
      while (*cc && *src != *cc) ++cc;
      return *cc ? src + 1 : 0;
    }

    // Match for members of char class.
    // Regex equivalent: /[axy]+/
    template <const char* char_class>
    const char* class_chars(const char* src) {
      const char* p = src;
      while (class_char<char_class>(p)) ++p;
      return p == src ? 0 : p;
    }

    // Match for members of char class.
    // Regex equivalent: /[^axy]/
    template <const char* neg_char_class>
    const char* neg_class_char(const char* src) {
      if (*src == 0) return 0;
      const char* cc = neg_char_class;
      while (*cc && *src != *cc) ++cc;
      return *cc ? 0 : src + 1;
    }

    // Match for members of char class.
    // Regex equivalent: /[^axy]+/
    template <const char* neg_char_class>
    const char* neg_class_chars(const char* src) {
      const char* p = src;
      while (neg_class_char<neg_char_class>(p)) ++p;
      return p == src ? 0 : p;
    }

    // Match all except the supplied one.
    // Regex equivalent: /[^x]/
    template <const char chr>
    const char* any_char_but(const char* src) {
      return (*src && *src != chr) ? src + 1 : 0;
    }

    // Succeeds if the matcher fails.
    // Aka. zero-width negative lookahead.
    // Regex equivalent: /(?!literal)/
    template <prelexer mx>
    const char* negate(const char* src) {
      return mx(src) ? 0 : src;
    }

    // Succeeds if the matcher succeeds.
    // Aka. zero-width positive lookahead.
    // Regex equivalent: /(?=literal)/
    // just hangs around until we need it
    template <prelexer mx>
    const char* lookahead(const char* src) {
      return mx(src) ? src : 0;
    }

    // Tries supplied matchers in order.
    // Succeeds if one of them succeeds.
    // Regex equivalent: /(?:FOO|BAR)/
    template <const prelexer mx>
    const char* alternatives(const char* src) {
      const char* rslt;
      if ((rslt = mx(src))) return rslt;
      return 0;
    }
    template <const prelexer mx1, const prelexer mx2, const prelexer... mxs>
    const char* alternatives(const char* src) {
      const char* rslt;
      if ((rslt = mx1(src))) return rslt;
      return alternatives<mx2, mxs...>(src);
    }

    // Tries supplied matchers in order.
    // Succeeds if all of them succeeds.
    // Regex equivalent: /(?:FOO)(?:BAR)/
    template <const prelexer mx1>
    const char* sequence(const char* src) {
      const char* rslt = src;
      if (!(rslt = mx1(rslt))) return 0;
      return rslt;
    }
    template <const prelexer mx1, const prelexer mx2, const prelexer... mxs>
    const char* sequence(const char* src) {
      const char* rslt = src;
      if (!(rslt = mx1(rslt))) return 0;
      return sequence<mx2, mxs...>(rslt);
    }


    // Match a pattern or not. Always succeeds.
    // Regex equivalent: /(?:literal)?/
    template <prelexer mx>
    const char* optional(const char* src) {
      const char* p = mx(src);
      return p ? p : src;
    }

    // Match zero or more of the patterns.
    // Regex equivalent: /(?:literal)*/
    template <prelexer mx>
    const char* zero_plus(const char* src) {
      const char* p = mx(src);
      while (p) src = p, p = mx(src);
      return src;
    }

    // Match one or more of the patterns.
    // Regex equivalent: /(?:literal)+/
    template <prelexer mx>
    const char* one_plus(const char* src) {
      const char* p = mx(src);
      if (!p) return 0;
      while (p) src = p, p = mx(src);
      return src;
    }

    // Match mx non-greedy until delimiter.
    // Other prelexers are greedy by default.
    // Regex equivalent: /(?:$mx)*?(?=$delim)\b/
    template <prelexer mx, prelexer delim>
    const char* non_greedy(const char* src) {
      while (!delim(src)) {
        const char* p = mx(src);
        if (p == src) return 0;
        if (p == 0) return 0;
        src = p;
      }
      return src;
    }

    //####################################
    // ADVANCED "REGEX" CONSTRUCTORS
    //####################################

    // Match with word boundary rule.
    // Regex equivalent: /(?:$mx)\b/i
    template <const char* str>
    const char* keyword(const char* src) {
      return sequence <
               insensitive < str >,
               word_boundary
             >(src);
    }

    // Match with word boundary rule.
    // Regex equivalent: /(?:$mx)\b/
    template <const char* str>
    const char* word(const char* src) {
      return sequence <
               exactly < str >,
               word_boundary
             >(src);
    }

    template <char chr>
    const char* loosely(const char* src) {
      return sequence <
               optional_spaces,
               exactly < chr >
             >(src);
    }
    template <const char* str>
    const char* loosely(const char* src) {
      return sequence <
               optional_spaces,
               exactly < str >
             >(src);
    }

  }
}

#endif

namespace Sass {
  // using namespace Lexer;
  namespace Prelexer {

    //####################################
    // KEYWORD "REGEX" MATCHERS
    //####################################

    // Match Sass boolean keywords.
    const char* kwd_true(const char* src);
    const char* kwd_false(const char* src);
    const char* kwd_only(const char* src);
    const char* kwd_and(const char* src);
    const char* kwd_or(const char* src);
    const char* kwd_not(const char* src);
    const char* kwd_eq(const char* src);
    const char* kwd_neq(const char* src);
    const char* kwd_gt(const char* src);
    const char* kwd_gte(const char* src);
    const char* kwd_lt(const char* src);
    const char* kwd_lte(const char* src);
    const char* kwd_using(const char* src);

    // Match standard control chars
    const char* kwd_at(const char* src);
    const char* kwd_dot(const char* src);
    const char* kwd_comma(const char* src);
    const char* kwd_colon(const char* src);
    const char* kwd_slash(const char* src);
    const char* kwd_star(const char* src);
    const char* kwd_plus(const char* src);
    const char* kwd_minus(const char* src);

    //####################################
    // SPECIAL "REGEX" CONSTRUCTS
    //####################################

    // Match a sequence of characters delimited by the supplied chars.
    template <char beg, char end, bool esc>
    const char* delimited_by(const char* src) {
      src = exactly<beg>(src);
      if (!src) return 0;
      const char* stop;
      while (true) {
        if (!*src) return 0;
        stop = exactly<end>(src);
        if (stop && (!esc || *(src - 1) != '\\')) return stop;
        src = stop ? stop : src + 1;
      }
    }

    // skip to delimiter (mx) inside given range
    // this will savely skip over all quoted strings
    // recursive skip stuff delimited by start/stop
    // first start/opener must be consumed already!
    template<prelexer start, prelexer stop>
    const char* skip_over_scopes(const char* src, const char* end) {

      size_t level = 0;
      bool in_squote = false;
      bool in_dquote = false;
      // bool in_braces = false;

      while (*src) {

        // check for abort condition
        if (end && src >= end) break;

        // has escaped sequence?
        if (*src == '\\') {
          ++ src; // skip this (and next)
        }
        else if (*src == '"') {
          in_dquote = ! in_dquote;
        }
        else if (*src == '\'') {
          in_squote = ! in_squote;
        }
        else if (in_dquote || in_squote) {
          // take everything literally
        }

        // find another opener inside?
        else if (const char* pos = start(src)) {
          ++ level; // increase counter
          src = pos - 1; // advance position
        }

        // look for the closer (maybe final, maybe not)
        else if (const char* final = stop(src)) {
          // only close one level?
          if (level > 0) -- level;
          // return position at end of stop
          // delimiter may be multiple chars
          else return final;
          // advance position
          src = final - 1;
        }

        // next
        ++ src;
      }

      return 0;
    }

    // skip to a skip delimited by parentheses
    // uses smart `skip_over_scopes` internally
    const char* parenthese_scope(const char* src);

    // skip to delimiter (mx) inside given range
    // this will savely skip over all quoted strings
    // recursive skip stuff delimited by start/stop
    // first start/opener must be consumed already!
    template<prelexer start, prelexer stop>
    const char* skip_over_scopes(const char* src) {
      return skip_over_scopes<start, stop>(src, 0);
    }

    // Match a sequence of characters delimited by the supplied chars.
    template <prelexer start, prelexer stop>
    const char* recursive_scopes(const char* src) {
      // parse opener
      src = start(src);
      // abort if not found
      if (!src) return 0;
      // parse the rest until final closer
      return skip_over_scopes<start, stop>(src);
    }

    // Match a sequence of characters delimited by the supplied strings.
    template <const char* beg, const char* end, bool esc>
    const char* delimited_by(const char* src) {
      src = exactly<beg>(src);
      if (!src) return 0;
      const char* stop;
      while (true) {
        if (!*src) return 0;
        stop = exactly<end>(src);
        if (stop && (!esc || *(src - 1) != '\\')) return stop;
        src = stop ? stop : src + 1;
      }
    }

    // Tries to match a certain number of times (between the supplied interval).
    template<prelexer mx, size_t lo, size_t hi>
    const char* between(const char* src) {
      for (size_t i = 0; i < lo; ++i) {
        src = mx(src);
        if (!src) return 0;
      }
      for (size_t i = lo; i <= hi; ++i) {
        const char* new_src = mx(src);
        if (!new_src) return src;
        src = new_src;
      }
      return src;
    }

    // equivalent of STRING_REGULAR_EXPRESSIONS
    const char* re_string_double_open(const char* src);
    const char* re_string_double_close(const char* src);
    const char* re_string_single_open(const char* src);
    const char* re_string_single_close(const char* src);
    const char* re_string_uri_open(const char* src);
    const char* re_string_uri_close(const char* src);

    // Match a line comment.
    const char* line_comment(const char* src);

    // Match a block comment.
    const char* block_comment(const char* src);
    // Match either.
    const char* comment(const char* src);
    // Match double- and single-quoted strings.
    const char* double_quoted_string(const char* src);
    const char* single_quoted_string(const char* src);
    const char* quoted_string(const char* src);
    // Match interpolants.
    const char* interpolant(const char* src);
    // Match number prefix ([\+\-]+)
    const char* number_prefix(const char* src);

    // Match zero plus white-space or line_comments
    const char* optional_css_whitespace(const char* src);
    const char* css_whitespace(const char* src);
    // Match optional_css_whitepace plus block_comments
    const char* optional_css_comments(const char* src);
    const char* css_comments(const char* src);

    // Match one backslash escaped char
    const char* escape_seq(const char* src);

    // Match CSS css variables.
    const char* custom_property_name(const char* src);
    // Match a CSS identifier.
    const char* identifier(const char* src);
    const char* identifier_alpha(const char* src);
    const char* identifier_alnum(const char* src);
    const char* strict_identifier(const char* src);
    const char* strict_identifier_alpha(const char* src);
    const char* strict_identifier_alnum(const char* src);
    // Match a CSS unit identifier.
    const char* one_unit(const char* src);
    const char* multiple_units(const char* src);
    const char* unit_identifier(const char* src);
    // const char* strict_identifier_alnums(const char* src);
    // Match reference selector.
    const char* re_reference_combinator(const char* src);
    const char* static_reference_combinator(const char* src);
    const char* schema_reference_combinator(const char* src);

    // Match interpolant schemas
    const char* identifier_schema(const char* src);
    const char* value_schema(const char* src);
    const char* sass_value(const char* src);
    // const char* filename(const char* src);
    // const char* filename_schema(const char* src);
    // const char* url_schema(const char* src);
    // const char* url_value(const char* src);
    const char* vendor_prefix(const char* src);

    const char* re_special_directive(const char* src);
    const char* re_prefixed_directive(const char* src);
    const char* re_almost_any_value_token(const char* src);

    // Match CSS '@' keywords.
    const char* at_keyword(const char* src);
    const char* kwd_import(const char* src);
    const char* kwd_at_root(const char* src);
    const char* kwd_with_directive(const char* src);
    const char* kwd_without_directive(const char* src);
    const char* kwd_media(const char* src);
    const char* kwd_supports_directive(const char* src);
    // const char* keyframes(const char* src);
    // const char* keyf(const char* src);
    const char* kwd_mixin(const char* src);
    const char* kwd_function(const char* src);
    const char* kwd_return_directive(const char* src);
    const char* kwd_include_directive(const char* src);
    const char* kwd_content_directive(const char* src);
    const char* kwd_charset_directive(const char* src);
    const char* kwd_extend(const char* src);

    const char* unicode_seq(const char* src);

    const char* kwd_if_directive(const char* src);
    const char* kwd_else_directive(const char* src);
    const char* elseif_directive(const char* src);

    const char* kwd_for_directive(const char* src);
    const char* kwd_from(const char* src);
    const char* kwd_to(const char* src);
    const char* kwd_through(const char* src);

    const char* kwd_each_directive(const char* src);
    const char* kwd_in(const char* src);

    const char* kwd_while_directive(const char* src);

    const char* re_nothing(const char* src);

    const char* re_special_fun(const char* src);

    const char* kwd_warn(const char* src);
    const char* kwd_err(const char* src);
    const char* kwd_dbg(const char* src);

    const char* kwd_null(const char* src);

    const char* re_selector_list(const char* src);
    const char* re_type_selector(const char* src);
    const char* re_static_expression(const char* src);

    // identifier that can start with hyphens
    const char* css_identifier(const char* src);
    const char* css_ip_identifier(const char* src);

    // Match CSS type selectors
    const char* namespace_schema(const char* src);
    const char* namespace_prefix(const char* src);
    const char* type_selector(const char* src);
    const char* hyphens_and_identifier(const char* src);
    const char* hyphens_and_name(const char* src);
    const char* universal(const char* src);
    // Match CSS id names.
    const char* id_name(const char* src);
    // Match CSS class names.
    const char* class_name(const char* src);
    // Attribute name in an attribute selector
    const char* attribute_name(const char* src);
    // Match placeholder selectors.
    const char* placeholder(const char* src);
    // Match CSS numeric constants.
    const char* op(const char* src);
    const char* sign(const char* src);
    const char* unsigned_number(const char* src);
    const char* number(const char* src);
    const char* coefficient(const char* src);
    const char* binomial(const char* src);
    const char* percentage(const char* src);
    const char* ampersand(const char* src);
    const char* dimension(const char* src);
    const char* hex(const char* src);
    const char* hexa(const char* src);
    const char* hex0(const char* src);
    // const char* rgb_prefix(const char* src);
    // Match CSS uri specifiers.
    const char* uri_prefix(const char* src);
    // Match CSS "!important" keyword.
    const char* kwd_important(const char* src);
    // Match CSS "!optional" keyword.
    const char* kwd_optional(const char* src);
    // Match Sass "!default" keyword.
    const char* default_flag(const char* src);
    const char* global_flag(const char* src);
    // Match CSS pseudo-class/element prefixes
    const char* pseudo_prefix(const char* src);
    // Match CSS function call openers.
    const char* re_functional(const char* src);
    const char* re_pseudo_selector(const char* src);
    const char* functional_schema(const char* src);
    const char* pseudo_not(const char* src);
    // Match CSS 'odd' and 'even' keywords for functional pseudo-classes.
    const char* even(const char* src);
    const char* odd(const char* src);
    // Match CSS attribute-matching operators.
    const char* exact_match(const char* src);
    const char* class_match(const char* src);
    const char* dash_match(const char* src);
    const char* prefix_match(const char* src);
    const char* suffix_match(const char* src);
    const char* substring_match(const char* src);
    // Match CSS combinators.
    // const char* adjacent_to(const char* src);
    // const char* precedes(const char* src);
    // const char* parent_of(const char* src);
    // const char* ancestor_of(const char* src);

    // Match SCSS variable names.
    const char* variable(const char* src);
    const char* calc_fn_call(const char* src);

    // IE stuff
    const char* ie_progid(const char* src);
    const char* ie_expression(const char* src);
    const char* ie_property(const char* src);
    const char* ie_keyword_arg(const char* src);
    const char* ie_keyword_arg_value(const char* src);
    const char* ie_keyword_arg_property(const char* src);

    // characters that terminate parsing of a list
    const char* list_terminator(const char* src);
    const char* space_list_terminator(const char* src);

    // match url()
    const char* H(const char* src);
    const char* W(const char* src);
    // `UNICODE` makes VS sad
    const char* UUNICODE(const char* src);
    const char* NONASCII(const char* src);
    const char* ESCAPE(const char* src);
    const char* real_uri(const char* src);
    const char* real_uri_suffix(const char* src);
    // const char* real_uri_prefix(const char* src);
    const char* real_uri_value(const char* src);

    // Path matching functions.
    // const char* folder(const char* src);
    // const char* folders(const char* src);


    const char* static_string(const char* src);
    const char* static_component(const char* src);
    const char* static_property(const char* src);
    const char* static_value(const char* src);

    const char* css_variable_value(const char* src);
    const char* css_variable_top_level_value(const char* src);

    // Utility functions for finding and counting characters in a string.
    template<char c>
    const char* find_first(const char* src) {
      while (*src && *src != c) ++src;
      return *src ? src : 0;
    }
    template<prelexer mx>
    const char* find_first(const char* src) {
      while (*src && !mx(src)) ++src;
      return *src ? src : 0;
    }
    template<prelexer mx>
    const char* find_first_in_interval(const char* beg, const char* end) {
      bool esc = false;
      while ((beg < end) && *beg) {
        if (esc) esc = false;
        else if (*beg == '\\') esc = true;
        else if (mx(beg)) return beg;
        ++beg;
      }
      return 0;
    }
    template<prelexer mx, prelexer skip>
    const char* find_first_in_interval(const char* beg, const char* end) {
      bool esc = false;
      while ((beg < end) && *beg) {
        if (esc) esc = false;
        else if (*beg == '\\') esc = true;
        else if (const char* pos = skip(beg)) beg = pos;
        else if (mx(beg)) return beg;
        ++beg;
      }
      return 0;
    }
    template <prelexer mx>
    unsigned int count_interval(const char* beg, const char* end) {
      unsigned int counter = 0;
      bool esc = false;
      while (beg < end && *beg) {
        const char* p;
        if (esc) {
          esc = false;
          ++beg;
        } else if (*beg == '\\') {
          esc = true;
          ++beg;
        } else if ((p = mx(beg))) {
          ++counter;
          beg = p;
        }
        else {
          ++beg;
        }
      }
      return counter;
    }

    template <size_t size, prelexer mx, prelexer pad>
    const char* padded_token(const char* src)
    {
      size_t got = 0;
      const char* pos = src;
      while (got < size) {
        if (!mx(pos)) break;
        ++ pos; ++ got;
      }
      while (got < size) {
        if (!pad(pos)) break;
        ++ pos; ++ got;
      }
      return got ? pos : 0;
    }

    template <size_t min, size_t max, prelexer mx>
    const char* minmax_range(const char* src)
    {
      size_t got = 0;
      const char* pos = src;
      while (got < max) {
        if (!mx(pos)) break;
        ++ pos; ++ got;
      }
      if (got < min) return 0;
      if (got > max) return 0;
      return pos;
    }

    template <char min, char max>
    const char* char_range(const char* src)
    {
      if (*src < min) return 0;
      if (*src > max) return 0;
      return src + 1;
    }

  }
}

#endif

#ifndef MAX_NESTING
// Note that this limit is not an exact science
// it depends on various factors, which some are
// not under our control (compile time or even OS
// dependent settings on the available stack size)
// It should fix most common segfault cases though.
#define MAX_NESTING 512
#endif

struct Lookahead {
  const char* found;
  const char* error;
  const char* position;
  bool parsable;
  bool has_interpolants;
  bool is_custom_property;
};

namespace Sass {

  class Parser : public ParserState {
  public:

    enum Scope { Root, Mixin, Function, Media, Control, Properties, Rules, AtRoot };

    Context& ctx;
    std::vector<Block_Obj> block_stack;
    std::vector<Scope> stack;
    Media_Block_Ptr last_media_block;
    const char* source;
    const char* position;
    const char* end;
    Position before_token;
    Position after_token;
    ParserState pstate;
    Backtraces traces;
    size_t indentation;
    size_t nestings;
    bool allow_parent;

    Token lexed;

    Parser(Context& ctx, const ParserState& pstate, Backtraces traces, bool allow_parent = true)
    : ParserState(pstate), ctx(ctx), block_stack(), stack(0), last_media_block(),
      source(0), position(0), end(0), before_token(pstate), after_token(pstate),
      pstate(pstate), traces(traces), indentation(0), nestings(0), allow_parent(allow_parent)
    {
      stack.push_back(Scope::Root);
    }

    // static Parser from_string(const std::string& src, Context& ctx, ParserState pstate = ParserState("[STRING]"));
    static Parser from_c_str(const char* src, Context& ctx, Backtraces, ParserState pstate = ParserState("[CSTRING]"), const char* source = nullptr, bool allow_parent = true);
    static Parser from_c_str(const char* beg, const char* end, Context& ctx, Backtraces, ParserState pstate = ParserState("[CSTRING]"), const char* source = nullptr, bool allow_parent = true);
    static Parser from_token(Token t, Context& ctx, Backtraces, ParserState pstate = ParserState("[TOKEN]"), const char* source = nullptr);
    // special static parsers to convert strings into certain selectors
    static Selector_List_Obj parse_selector(const char* src, Context& ctx, Backtraces, ParserState pstate = ParserState("[SELECTOR]"), const char* source = nullptr, bool allow_parent = true);

#ifdef __clang__

    // lex and peak uses the template parameter to branch on the action, which
    // triggers clangs tautological comparison on the single-comparison
    // branches. This is not a bug, just a merging of behaviour into
    // one function

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wtautological-compare"

#endif


    // skip current token and next whitespace
    // moves ParserState right before next token
    void advanceToNextToken();

    bool peek_newline(const char* start = 0);

    // skip over spaces, tabs and line comments
    template <Prelexer::prelexer mx>
    const char* sneak(const char* start = 0)
    {
      using namespace Prelexer;

      // maybe use optional start position from arguments?
      const char* it_position = start ? start : position;

      // skip white-space?
      if (mx == spaces ||
          mx == no_spaces ||
          mx == css_comments ||
          mx == css_whitespace ||
          mx == optional_spaces ||
          mx == optional_css_comments ||
          mx == optional_css_whitespace
      ) {
        return it_position;
      }

      // skip over spaces, tabs and sass line comments
      const char* pos = optional_css_whitespace(it_position);
      // always return a valid position
      return pos ? pos : it_position;

    }

    // match will not skip over space, tabs and line comment
    // return the position where the lexer match will occur
    template <Prelexer::prelexer mx>
    const char* match(const char* start = 0)
    {
      // match the given prelexer
      return mx(position);
    }

    // peek will only skip over space, tabs and line comment
    // return the position where the lexer match will occur
    template <Prelexer::prelexer mx>
    const char* peek(const char* start = 0)
    {

      // sneak up to the actual token we want to lex
      // this should skip over white-space if desired
      const char* it_before_token = sneak < mx >(start);

      // match the given prelexer
      const char* match = mx(it_before_token);

      // check if match is in valid range
      return match <= end ? match : 0;

    }

    // white-space handling is built into the lexer
    // this way you do not need to parse it yourself
    // some matchers don't accept certain white-space
    // we do not support start arg, since we manipulate
    // sourcemap offset and we modify the position pointer!
    // lex will only skip over space, tabs and line comment
    template <Prelexer::prelexer mx>
    const char* lex(bool lazy = true, bool force = false)
    {

      if (*position == 0) return 0;

      // position considered before lexed token
      // we can skip whitespace or comments for
      // lazy developers (but we need control)
      const char* it_before_token = position;

      // sneak up to the actual token we want to lex
      // this should skip over white-space if desired
      if (lazy) it_before_token = sneak < mx >(position);

      // now call matcher to get position after token
      const char* it_after_token = mx(it_before_token);

      // check if match is in valid range
      if (it_after_token > end) return 0;

      // maybe we want to update the parser state anyway?
      if (force == false) {
        // assertion that we got a valid match
        if (it_after_token == 0) return 0;
        // assertion that we actually lexed something
        if (it_after_token == it_before_token) return 0;
      }

      // create new lexed token object (holds the parse results)
      lexed = Token(position, it_before_token, it_after_token);

      // advance position (add whitespace before current token)
      before_token = after_token.add(position, it_before_token);

      // update after_token position for current token
      after_token.add(it_before_token, it_after_token);

      // ToDo: could probably do this incremetal on original object (API wants offset?)
      pstate = ParserState(path, source, lexed, before_token, after_token - before_token);

      // advance internal char iterator
      return position = it_after_token;

    }

    // lex_css skips over space, tabs, line and block comment
    // all block comments will be consumed and thrown away
    // source-map position will point to token after the comment
    template <Prelexer::prelexer mx>
    const char* lex_css()
    {
      // copy old token
      Token prev = lexed;
      // store previous pointer
      const char* oldpos = position;
      Position bt = before_token;
      Position at = after_token;
      ParserState op = pstate;
      // throw away comments
      // update srcmap position
      lex < Prelexer::css_comments >();
      // now lex a new token
      const char* pos = lex< mx >();
      // maybe restore prev state
      if (pos == 0) {
        pstate = op;
        lexed = prev;
        position = oldpos;
        after_token = at;
        before_token = bt;
      }
      // return match
      return pos;
    }

    // all block comments will be skipped and thrown away
    template <Prelexer::prelexer mx>
    const char* peek_css(const char* start = 0)
    {
      // now peek a token (skip comments first)
      return peek< mx >(peek < Prelexer::css_comments >(start));
    }

#ifdef __clang__

#pragma clang diagnostic pop

#endif

    void error(std::string msg);
    void error(std::string msg, Position pos);
    // generate message with given and expected sample
    // text before and in the middle are configurable
    void css_error(const std::string& msg,
                   const std::string& prefix = " after ",
                   const std::string& middle = ", was: ",
                   const bool trim = true);
    void read_bom();

    Block_Obj parse();
    Import_Obj parse_import();
    Definition_Obj parse_definition(Definition::Type which_type);
    Parameters_Obj parse_parameters();
    Parameter_Obj parse_parameter();
    Mixin_Call_Obj parse_include_directive();
    Arguments_Obj parse_arguments();
    Argument_Obj parse_argument();
    Assignment_Obj parse_assignment();
    Ruleset_Obj parse_ruleset(Lookahead lookahead);
    Selector_List_Obj parse_selector_list(bool chroot);
    Complex_Selector_Obj parse_complex_selector(bool chroot);
    Selector_Schema_Obj parse_selector_schema(const char* end_of_selector, bool chroot);
    Compound_Selector_Obj parse_compound_selector();
    Simple_Selector_Obj parse_simple_selector();
    Wrapped_Selector_Obj parse_negated_selector();
    Simple_Selector_Obj parse_pseudo_selector();
    Attribute_Selector_Obj parse_attribute_selector();
    Block_Obj parse_block(bool is_root = false);
    Block_Obj parse_css_block(bool is_root = false);
    bool parse_block_nodes(bool is_root = false);
    bool parse_block_node(bool is_root = false);

    bool parse_number_prefix();
    Declaration_Obj parse_declaration();
    Expression_Obj parse_map();
    Expression_Obj parse_bracket_list();
    Expression_Obj parse_list(bool delayed = false);
    Expression_Obj parse_comma_list(bool delayed = false);
    Expression_Obj parse_space_list();
    Expression_Obj parse_disjunction();
    Expression_Obj parse_conjunction();
    Expression_Obj parse_relation();
    Expression_Obj parse_expression();
    Expression_Obj parse_operators();
    Expression_Obj parse_factor();
    Expression_Obj parse_value();
    Function_Call_Obj parse_calc_function();
    Function_Call_Obj parse_function_call();
    Function_Call_Obj parse_function_call_schema();
    String_Obj parse_url_function_string();
    String_Obj parse_url_function_argument();
    String_Obj parse_interpolated_chunk(Token, bool constant = false, bool css = true);
    String_Obj parse_string();
    Value_Obj parse_static_value();
    String_Schema_Obj parse_css_variable_value(bool top_level = true);
    String_Schema_Obj parse_css_variable_value_token(bool top_level = true);
    String_Obj parse_ie_property();
    String_Obj parse_ie_keyword_arg();
    String_Schema_Obj parse_value_schema(const char* stop);
    String_Obj parse_identifier_schema();
    If_Obj parse_if_directive(bool else_if = false);
    For_Obj parse_for_directive();
    Each_Obj parse_each_directive();
    While_Obj parse_while_directive();
    Return_Obj parse_return_directive();
    Content_Obj parse_content_directive();
    void parse_charset_directive();
    Media_Block_Obj parse_media_block();
    List_Obj parse_media_queries();
    Media_Query_Obj parse_media_query();
    Media_Query_Expression_Obj parse_media_expression();
    Supports_Block_Obj parse_supports_directive();
    Supports_Condition_Obj parse_supports_condition();
    Supports_Condition_Obj parse_supports_negation();
    Supports_Condition_Obj parse_supports_operator();
    Supports_Condition_Obj parse_supports_interpolation();
    Supports_Condition_Obj parse_supports_declaration();
    Supports_Condition_Obj parse_supports_condition_in_parens();
    At_Root_Block_Obj parse_at_root_block();
    At_Root_Query_Obj parse_at_root_query();
    String_Schema_Obj parse_almost_any_value();
    Directive_Obj parse_special_directive();
    Directive_Obj parse_prefixed_directive();
    Directive_Obj parse_directive();
    Warning_Obj parse_warning();
    Error_Obj parse_error();
    Debug_Obj parse_debug();

    Value_Ptr color_or_string(const std::string& lexed) const;

    // be more like ruby sass
    Expression_Obj lex_almost_any_value_token();
    Expression_Obj lex_almost_any_value_chars();
    Expression_Obj lex_interp_string();
    Expression_Obj lex_interp_uri();
    Expression_Obj lex_interpolation();

    // these will throw errors
    Token lex_variable();
    Token lex_identifier();

    void parse_block_comments();

    Lookahead lookahead_for_value(const char* start = 0);
    Lookahead lookahead_for_selector(const char* start = 0);
    Lookahead lookahead_for_include(const char* start = 0);

    Expression_Obj fold_operands(Expression_Obj base, std::vector<Expression_Obj>& operands, Operand op);
    Expression_Obj fold_operands(Expression_Obj base, std::vector<Expression_Obj>& operands, std::vector<Operand>& ops, size_t i = 0);

    void throw_syntax_error(std::string message, size_t ln = 0);
    void throw_read_error(std::string message, size_t ln = 0);


    template <Prelexer::prelexer open, Prelexer::prelexer close>
    Expression_Obj lex_interp()
    {
      if (lex < open >(false)) {
        String_Schema_Obj schema = SASS_MEMORY_NEW(String_Schema, pstate);
        // std::cerr << "LEX [[" << std::string(lexed) << "]]\n";
        schema->append(SASS_MEMORY_NEW(String_Constant, pstate, lexed));
        if (position[0] == '#' && position[1] == '{') {
          Expression_Obj itpl = lex_interpolation();
          if (!itpl.isNull()) schema->append(itpl);
          while (lex < close >(false)) {
            // std::cerr << "LEX [[" << std::string(lexed) << "]]\n";
            schema->append(SASS_MEMORY_NEW(String_Constant, pstate, lexed));
            if (position[0] == '#' && position[1] == '{') {
              Expression_Obj itpl = lex_interpolation();
              if (!itpl.isNull()) schema->append(itpl);
            } else {
              return schema;
            }
          }
        } else {
          return SASS_MEMORY_NEW(String_Constant, pstate, lexed);
        }
      }
      return {};
    }

  public:
    static Number_Ptr lexed_number(const ParserState& pstate, const std::string& parsed);
    static Number_Ptr lexed_dimension(const ParserState& pstate, const std::string& parsed);
    static Number_Ptr lexed_percentage(const ParserState& pstate, const std::string& parsed);
    static Value_Ptr lexed_hex_color(const ParserState& pstate, const std::string& parsed);
  private:
    Number_Ptr lexed_number(const std::string& parsed) { return lexed_number(pstate, parsed); };
    Number_Ptr lexed_dimension(const std::string& parsed) { return lexed_dimension(pstate, parsed); };
    Number_Ptr lexed_percentage(const std::string& parsed) { return lexed_percentage(pstate, parsed); };
    Value_Ptr lexed_hex_color(const std::string& parsed) { return lexed_hex_color(pstate, parsed); };

    static const char* re_attr_sensitive_close(const char* src);
    static const char* re_attr_insensitive_close(const char* src);

  };

  size_t check_bom_chars(const char* src, const char *end, const unsigned char* bom, size_t len);
}

#endif
#ifndef SASS_SASS_UTIL_H
#define SASS_SASS_UTIL_H


namespace Sass {




  /*
   This is for ports of functions in the Sass:Util module.
   */


  /*
    # Return a Node collection of all possible paths through the given Node collection of Node collections.
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
  */
  Node paths(const Node& arrs);


  /*
  This class is a default implementation of a Node comparator that can be passed to the lcs function below.
  It uses operator== for equality comparision. It then returns one if the Nodes are equal.
  */
  class DefaultLcsComparator {
  public:
    bool operator()(const Node& one, const Node& two, Node& out) const {
      // TODO: Is this the correct C++ interpretation?
      // block ||= proc {|a, b| a == b && a}
      if (one == two) {
        out = one;
        return true;
      }

      return false;
    }
  };


  typedef std::vector<std::vector<int> > LCSTable;


  /*
  This is the equivalent of ruby's Sass::Util.lcs_backtrace.

  # Computes a single longest common subsequence for arrays x and y.
  # Algorithm from http://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Reading_out_an_LCS
  */
  template<typename ComparatorType>
  Node lcs_backtrace(const LCSTable& c, const Node& x, const Node& y, int i, int j, const ComparatorType& comparator) {
    DEBUG_PRINTLN(LCS, "LCSBACK: X=" << x << " Y=" << y << " I=" << i << " J=" << j)

    if (i == 0 || j == 0) {
      DEBUG_PRINTLN(LCS, "RETURNING EMPTY")
      return Node::createCollection();
    }

    NodeDeque& xChildren = *(x.collection());
    NodeDeque& yChildren = *(y.collection());

    Node compareOut = Node::createNil();
    if (comparator(xChildren[i], yChildren[j], compareOut)) {
      DEBUG_PRINTLN(LCS, "RETURNING AFTER ELEM COMPARE")
      Node result = lcs_backtrace(c, x, y, i - 1, j - 1, comparator);
      result.collection()->push_back(compareOut);
      return result;
    }

    if (c[i][j - 1] > c[i - 1][j]) {
      DEBUG_PRINTLN(LCS, "RETURNING AFTER TABLE COMPARE")
      return lcs_backtrace(c, x, y, i, j - 1, comparator);
    }

    DEBUG_PRINTLN(LCS, "FINAL RETURN")
    return lcs_backtrace(c, x, y, i - 1, j, comparator);
  }


  /*
  This is the equivalent of ruby's Sass::Util.lcs_table.

  # Calculates the memoization table for the Least Common Subsequence algorithm.
  # Algorithm from http://en.wikipedia.org/wiki/Longest_common_subsequence_problem#Computing_the_length_of_the_LCS
  */
  template<typename ComparatorType>
  void lcs_table(const Node& x, const Node& y, const ComparatorType& comparator, LCSTable& out) {
    DEBUG_PRINTLN(LCS, "LCSTABLE: X=" << x << " Y=" << y)

    NodeDeque& xChildren = *(x.collection());
    NodeDeque& yChildren = *(y.collection());

    LCSTable c(xChildren.size(), std::vector<int>(yChildren.size()));

    // These shouldn't be necessary since the vector will be initialized to 0 already.
    // x.size.times {|i| c[i][0] = 0}
    // y.size.times {|j| c[0][j] = 0}

    for (size_t i = 1; i < xChildren.size(); i++) {
      for (size_t j = 1; j < yChildren.size(); j++) {
        Node compareOut = Node::createNil();

        if (comparator(xChildren[i], yChildren[j], compareOut)) {
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
  template<typename ComparatorType>
  Node lcs(Node& x, Node& y, const ComparatorType& comparator) {
    DEBUG_PRINTLN(LCS, "LCS: X=" << x << " Y=" << y)

    Node newX = Node::createCollection();
    newX.collection()->push_back(Node::createNil());
    newX.plus(x);

    Node newY = Node::createCollection();
    newY.collection()->push_back(Node::createNil());
    newY.plus(y);

    LCSTable table;
    lcs_table(newX, newY, comparator, table);

    return lcs_backtrace(table, newX, newY, static_cast<int>(newX.collection()->size()) - 1, static_cast<int>(newY.collection()->size()) - 1, comparator);
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
  */
  Node flatten(Node& arr, int n = -1);


  /*
  This is the equivalent of ruby's Sass::Util.group_by_to_a.

  # Performs the equivalent of `enum.group_by.to_a`, but with a guaranteed
  # order. Unlike [#hash_to_a], the resulting order isn't sorted key order;
  # instead, it's the same order as `#group_by` has under Ruby 1.9 (key
  # appearance order).
  #
  # @param enum [Enumerable]
  # @return [Array<[Object, Array]>] An array of pairs.

  TODO: update @param and @return once I know what those are.

  The following is the modified version of the ruby code that was more portable to C++. You
  should be able to drop it into ruby 3.2.19 and get the same results from ruby sass.

    def group_by_to_a(enum, &block)
      order = {}

      arr = []

      grouped = {}

      for e in enum do
        key = block[e]
        unless order.include?(key)
          order[key] = order.size
        end

        if not grouped.has_key?(key) then
          grouped[key] = [e]
        else
          grouped[key].push(e)
        end
      end

      grouped.each do |key, vals|
        arr[order[key]] = [key, vals]
      end

      arr
    end

  */
  template<typename EnumType, typename KeyType, typename KeyFunctorType>
  void group_by_to_a(std::vector<EnumType>& enumeration, KeyFunctorType& keyFunc, std::vector<std::pair<KeyType, std::vector<EnumType> > >& arr /*out*/) {

    std::map<unsigned int, KeyType> order;

    std::map<size_t, std::vector<EnumType> > grouped;

    for (typename std::vector<EnumType>::iterator enumIter = enumeration.begin(), enumIterEnd = enumeration.end(); enumIter != enumIterEnd; enumIter++) {
      EnumType& e = *enumIter;

      KeyType key = keyFunc(e);

      if (grouped.find(key->hash()) == grouped.end()) {
        order.insert(std::make_pair((unsigned int)order.size(), key));

        std::vector<EnumType> newCollection;
        newCollection.push_back(e);
        grouped.insert(std::make_pair(key->hash(), newCollection));
      } else {
        std::vector<EnumType>& collection = grouped.at(key->hash());
        collection.push_back(e);
      }
    }

    for (unsigned int index = 0; index < order.size(); index++) {
      KeyType& key = order.at(index);
      std::vector<EnumType>& values = grouped.at(key->hash());

      std::pair<KeyType, std::vector<EnumType> > grouping = std::make_pair(key, values);

      arr.push_back(grouping);
    }
  }


}

#endif
#ifndef SASS_CHECK_NESTING_H
#define SASS_CHECK_NESTING_H


namespace Sass {

  class CheckNesting : public Operation_CRTP<Statement_Ptr, CheckNesting> {

    std::vector<Statement_Ptr>  parents;
    Backtraces                  traces;
    Statement_Ptr               parent;
    Definition_Ptr              current_mixin_definition;

    Statement_Ptr before(Statement_Ptr);
    Statement_Ptr visit_children(Statement_Ptr);

  public:
    CheckNesting();
    ~CheckNesting() { }

    Statement_Ptr operator()(Block_Ptr);
    Statement_Ptr operator()(Definition_Ptr);
    Statement_Ptr operator()(If_Ptr);

    template <typename U>
    Statement_Ptr fallback(U x) {
      Statement_Ptr s = Cast<Statement>(x);
      if (s && this->should_visit(s)) {
        Block_Ptr b1 = Cast<Block>(s);
        Has_Block_Ptr b2 = Cast<Has_Block>(s);
        if (b1 || b2) return visit_children(s);
      }
      return s;
    }

  private:
    void invalid_content_parent(Statement_Ptr, AST_Node_Ptr);
    void invalid_charset_parent(Statement_Ptr, AST_Node_Ptr);
    void invalid_extend_parent(Statement_Ptr, AST_Node_Ptr);
    // void invalid_import_parent(Statement_Ptr);
    void invalid_mixin_definition_parent(Statement_Ptr, AST_Node_Ptr);
    void invalid_function_parent(Statement_Ptr, AST_Node_Ptr);

    void invalid_function_child(Statement_Ptr);
    void invalid_prop_child(Statement_Ptr);
    void invalid_prop_parent(Statement_Ptr, AST_Node_Ptr);
    void invalid_return_parent(Statement_Ptr, AST_Node_Ptr);
    void invalid_value_child(AST_Node_Ptr);

    bool is_transparent_parent(Statement_Ptr, Statement_Ptr);

    bool should_visit(Statement_Ptr);

    bool is_charset(Statement_Ptr);
    bool is_mixin(Statement_Ptr);
    bool is_function(Statement_Ptr);
    bool is_root_node(Statement_Ptr);
    bool is_at_root_node(Statement_Ptr);
    bool is_directive_node(Statement_Ptr);
  };

}

#endif
#ifndef SASS_FN_SELECTORS_H
#define SASS_FN_SELECTORS_H


namespace Sass {

  namespace Functions {

    #define ARGSEL(argname) get_arg_sel(argname, env, sig, pstate, traces, ctx)
    #define ARGSELS(argname) get_arg_sels(argname, env, sig, pstate, traces, ctx)

    BUILT_IN(selector_nest);
    BUILT_IN(selector_append);
    BUILT_IN(selector_extend);
    BUILT_IN(selector_replace);
    BUILT_IN(selector_unify);
    BUILT_IN(is_superselector);
    BUILT_IN(simple_selectors);
    BUILT_IN(selector_parse);

    extern Signature selector_nest_sig;
    extern Signature selector_append_sig;
    extern Signature selector_extend_sig;
    extern Signature selector_replace_sig;
    extern Signature selector_unify_sig;
    extern Signature is_superselector_sig;
    extern Signature simple_selectors_sig;
    extern Signature selector_parse_sig;

  }

}

#endif
#ifndef SASS_CSSIZE_H
#define SASS_CSSIZE_H


namespace Sass {

  struct Backtrace;

  class Cssize : public Operation_CRTP<Statement_Ptr, Cssize> {

    Context&                    ctx;
    Backtraces&                 traces;
    BlockStack      block_stack;
    std::vector<Statement_Ptr>  p_stack;

  public:
    Cssize(Context&);
    ~Cssize() { }

    Selector_List_Ptr selector();

    Block_Ptr operator()(Block_Ptr);
    Statement_Ptr operator()(Ruleset_Ptr);
    // Statement_Ptr operator()(Bubble_Ptr);
    Statement_Ptr operator()(Media_Block_Ptr);
    Statement_Ptr operator()(Supports_Block_Ptr);
    Statement_Ptr operator()(At_Root_Block_Ptr);
    Statement_Ptr operator()(Directive_Ptr);
    Statement_Ptr operator()(Keyframe_Rule_Ptr);
    Statement_Ptr operator()(Trace_Ptr);
    Statement_Ptr operator()(Declaration_Ptr);
    // Statement_Ptr operator()(Assignment_Ptr);
    // Statement_Ptr operator()(Import_Ptr);
    // Statement_Ptr operator()(Import_Stub_Ptr);
    // Statement_Ptr operator()(Warning_Ptr);
    // Statement_Ptr operator()(Error_Ptr);
    // Statement_Ptr operator()(Comment_Ptr);
    // Statement_Ptr operator()(If_Ptr);
    // Statement_Ptr operator()(For_Ptr);
    // Statement_Ptr operator()(Each_Ptr);
    // Statement_Ptr operator()(While_Ptr);
    // Statement_Ptr operator()(Return_Ptr);
    // Statement_Ptr operator()(Extension_Ptr);
    // Statement_Ptr operator()(Definition_Ptr);
    // Statement_Ptr operator()(Mixin_Call_Ptr);
    // Statement_Ptr operator()(Content_Ptr);
    Statement_Ptr operator()(Null_Ptr);

    Statement_Ptr parent();
    std::vector<std::pair<bool, Block_Obj>> slice_by_bubble(Block_Ptr);
    Statement_Ptr bubble(Directive_Ptr);
    Statement_Ptr bubble(At_Root_Block_Ptr);
    Statement_Ptr bubble(Media_Block_Ptr);
    Statement_Ptr bubble(Supports_Block_Ptr);

    Block_Ptr debubble(Block_Ptr children, Statement_Ptr parent = 0);
    Block_Ptr flatten(Block_Ptr);
    bool bubblable(Statement_Ptr);

    List_Ptr merge_media_queries(Media_Block_Ptr, Media_Block_Ptr);
    Media_Query_Ptr merge_media_query(Media_Query_Ptr, Media_Query_Ptr);

    // generic fallback
    template <typename U>
    Statement_Ptr fallback(U x)
    { return Cast<Statement>(x); }

    void append_block(Block_Ptr, Block_Ptr);
  };

}

#endif
#ifndef SASS_PATHS_H
#define SASS_PATHS_H



template<typename T>
std::string vector_to_string(std::vector<T> v)
{
  std::stringstream buffer;
  buffer << "[";

  if (!v.empty())
  {  buffer << v[0]; }
  else
  { buffer << "]"; }

  if (v.size() == 1)
  { buffer << "]"; }
  else
  {
    for (size_t i = 1, S = v.size(); i < S; ++i) buffer << ", " << v[i];
    buffer << "]";
  }

  return buffer.str();
}

namespace Sass {


  template<typename T>
  std::vector<std::vector<T> > paths(std::vector<std::vector<T> > strata, size_t from_end = 0)
  {
    if (strata.empty()) {
      return std::vector<std::vector<T> >();
    }

    size_t end = strata.size() - from_end;
    if (end <= 1) {
      std::vector<std::vector<T> > starting_points;
      starting_points.reserve(strata[0].size());
      for (size_t i = 0, S = strata[0].size(); i < S; ++i) {
        std::vector<T> starting_point;
        starting_point.push_back(strata[0][i]);
        starting_points.push_back(starting_point);
      }
      return starting_points;
    }

    std::vector<std::vector<T> > up_to_here = paths(strata, from_end + 1);
    std::vector<T>          here       = strata[end-1];

    std::vector<std::vector<T> > branches;
    branches.reserve(up_to_here.size() * here.size());
    for (size_t i = 0, S1 = up_to_here.size(); i < S1; ++i) {
      for (size_t j = 0, S2 = here.size(); j < S2; ++j) {
        std::vector<T> branch = up_to_here[i];
        branch.push_back(here[j]);
        branches.push_back(branch);
      }
    }

    return branches;
  }

}

#endif
#ifndef SASS_FN_LISTS_H
#define SASS_FN_LISTS_H


namespace Sass {

  namespace Functions {

    extern Signature length_sig;
    extern Signature nth_sig;
    extern Signature index_sig;
    extern Signature join_sig;
    extern Signature append_sig;
    extern Signature zip_sig;
    extern Signature list_separator_sig;
    extern Signature is_bracketed_sig;
    extern Signature keywords_sig;

    BUILT_IN(length);
    BUILT_IN(nth);
    BUILT_IN(index);
    BUILT_IN(join);
    BUILT_IN(append);
    BUILT_IN(zip);
    BUILT_IN(list_separator);
    BUILT_IN(is_bracketed);
    BUILT_IN(keywords);

  }

}

#endif
#ifndef SASS_DEBUGGER_H
#define SASS_DEBUGGER_H


using namespace Sass;

inline void debug_ast(AST_Node_Ptr node, std::string ind = "", Env* env = 0);

inline void debug_ast(const AST_Node* node, std::string ind = "", Env* env = 0) {
  debug_ast(const_cast<AST_Node*>(node), ind, env);
}

inline void debug_sources_set(ComplexSelectorSet& set, std::string ind = "")
{
  if (ind == "") std::cerr << "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  for(auto const &pair : set) {
    debug_ast(pair, ind + "");
    // debug_ast(set[pair], ind + "first: ");
  }
  if (ind == "") std::cerr << "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
}

inline std::string str_replace(std::string str, const std::string& oldStr, const std::string& newStr)
{
  size_t pos = 0;
  while((pos = str.find(oldStr, pos)) != std::string::npos)
  {
     str.replace(pos, oldStr.length(), newStr);
     pos += newStr.length();
  }
  return str;
}

inline std::string prettyprint(const std::string& str) {
  std::string clean = str_replace(str, "\n", "\\n");
  clean = str_replace(clean, "	", "\\t");
  clean = str_replace(clean, "\r", "\\r");
  return clean;
}

inline std::string longToHex(long long t) {
  std::stringstream is;
  is << std::hex << t;
  return is.str();
}

inline std::string pstate_source_position(AST_Node_Ptr node)
{
  std::stringstream str;
  Position start(node->pstate());
  Position end(start + node->pstate().offset);
  str << (start.file == std::string::npos ? -1 : start.file)
    << "@[" << start.line << ":" << start.column << "]"
    << "-[" << end.line << ":" << end.column << "]";
#ifdef DEBUG_SHARED_PTR
      str << "x" << node->getRefCount() << ""
      << " " << node->getDbgFile()
      << "@" << node->getDbgLine();
#endif
  return str.str();
}

inline void debug_ast(AST_Node_Ptr node, std::string ind, Env* env)
{
  if (node == 0) return;
  if (ind == "") std::cerr << "####################################################################\n";
  if (Cast<Bubble>(node)) {
    Bubble_Ptr bubble = Cast<Bubble>(node);
    std::cerr << ind << "Bubble " << bubble;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << bubble->tabs();
    std::cerr << std::endl;
    debug_ast(bubble->node(), ind + " ", env);
  } else if (Cast<Trace>(node)) {
    Trace_Ptr trace = Cast<Trace>(node);
    std::cerr << ind << "Trace " << trace;
    std::cerr << " (" << pstate_source_position(node) << ")"
    << " [name:" << trace->name() << ", type: " << trace->type() << "]"
    << std::endl;
    debug_ast(trace->block(), ind + " ", env);
  } else if (Cast<At_Root_Block>(node)) {
    At_Root_Block_Ptr root_block = Cast<At_Root_Block>(node);
    std::cerr << ind << "At_Root_Block " << root_block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << root_block->tabs();
    std::cerr << std::endl;
    debug_ast(root_block->expression(), ind + ":", env);
    debug_ast(root_block->block(), ind + " ", env);
  } else if (Cast<Selector_List>(node)) {
    Selector_List_Ptr selector = Cast<Selector_List>(node);
    std::cerr << ind << "Selector_List " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " [@media:" << selector->media_block() << "]";
    std::cerr << (selector->is_invisible() ? " [INVISIBLE]": " -");
    std::cerr << (selector->has_placeholder() ? " [PLACEHOLDER]": " -");
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << std::endl;
    debug_ast(selector->schema(), ind + "#{} ");

    for(const Complex_Selector_Obj& i : selector->elements()) { debug_ast(i, ind + " ", env); }

//  } else if (Cast<Expression>(node)) {
//    Expression_Ptr expression = Cast<Expression>(node);
//    std::cerr << ind << "Expression " << expression << " " << expression->concrete_type() << std::endl;

  } else if (Cast<Parent_Reference>(node)) {
    Parent_Reference_Ptr selector = Cast<Parent_Reference>(node);
    std::cerr << ind << "Parent_Reference " << selector;
//    if (selector->not_selector()) cerr << " [in_declaration]";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <" << prettyprint(selector->pstate().token.ws_before()) << ">" << std::endl;
//    debug_ast(selector->selector(), ind + "->", env);

  } else if (Cast<Parent_Selector>(node)) {
    Parent_Selector_Ptr selector = Cast<Parent_Selector>(node);
    std::cerr << ind << "Parent_Selector " << selector;
//    if (selector->not_selector()) cerr << " [in_declaration]";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " [" << (selector->real() ? "REAL" : "FAKE") << "]";
    std::cerr << " <" << prettyprint(selector->pstate().token.ws_before()) << ">" << std::endl;
//    debug_ast(selector->selector(), ind + "->", env);

  } else if (Cast<Complex_Selector>(node)) {
    Complex_Selector_Ptr selector = Cast<Complex_Selector>(node);
    std::cerr << ind << "Complex_Selector " << selector
      << " (" << pstate_source_position(node) << ")"
      << " <" << selector->hash() << ">"
      << " [length:" << longToHex(selector->length()) << "]"
      << " [weight:" << longToHex(selector->specificity()) << "]"
      << " [@media:" << selector->media_block() << "]"
      << (selector->is_invisible() ? " [INVISIBLE]": " -")
      << (selector->has_placeholder() ? " [PLACEHOLDER]": " -")
      << (selector->is_optional() ? " [is_optional]": " -")
      << (selector->has_parent_ref() ? " [has parent]": " -")
      << (selector->has_line_feed() ? " [line-feed]": " -")
      << (selector->has_line_break() ? " [line-break]": " -")
      << " -- ";
      std::string del;
      switch (selector->combinator()) {
        case Complex_Selector::PARENT_OF:   del = ">"; break;
        case Complex_Selector::PRECEDES:    del = "~"; break;
        case Complex_Selector::ADJACENT_TO: del = "+"; break;
        case Complex_Selector::ANCESTOR_OF: del = " "; break;
        case Complex_Selector::REFERENCE:   del = "//"; break;
      }
      // if (del = "/") del += selector->reference()->perform(&to_string) + "/";
    std::cerr << " <" << prettyprint(selector->pstate().token.ws_before()) << ">" << std::endl;
    debug_ast(selector->head(), ind + " " /* + "[" + del + "]" */, env);
    if (selector->tail()) {
      debug_ast(selector->tail(), ind + "{" + del + "}", env);
    } else if(del != " ") {
      std::cerr << ind << " |" << del << "| {trailing op}" << std::endl;
    }
    ComplexSelectorSet set = selector->sources();
    // debug_sources_set(set, ind + "  @--> ");
  } else if (Cast<Compound_Selector>(node)) {
    Compound_Selector_Ptr selector = Cast<Compound_Selector>(node);
    std::cerr << ind << "Compound_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " [weight:" << longToHex(selector->specificity()) << "]";
    std::cerr << " [@media:" << selector->media_block() << "]";
    std::cerr << (selector->extended() ? " [extended]": " -");
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << " <" << prettyprint(selector->pstate().token.ws_before()) << ">" << std::endl;
    for(const Simple_Selector_Obj& i : selector->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Wrapped_Selector>(node)) {
    Wrapped_Selector_Ptr selector = Cast<Wrapped_Selector>(node);
    std::cerr << ind << "Wrapped_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <<" << selector->ns_name() << ">>";
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << std::endl;
    debug_ast(selector->selector(), ind + " () ", env);
  } else if (Cast<Pseudo_Selector>(node)) {
    Pseudo_Selector_Ptr selector = Cast<Pseudo_Selector>(node);
    std::cerr << ind << "Pseudo_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <<" << selector->ns_name() << ">>";
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << std::endl;
    debug_ast(selector->expression(), ind + " <= ", env);
  } else if (Cast<Attribute_Selector>(node)) {
    Attribute_Selector_Ptr selector = Cast<Attribute_Selector>(node);
    std::cerr << ind << "Attribute_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <<" << selector->ns_name() << ">>";
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << std::endl;
    debug_ast(selector->value(), ind + "[" + selector->matcher() + "] ", env);
  } else if (Cast<Class_Selector>(node)) {
    Class_Selector_Ptr selector = Cast<Class_Selector>(node);
    std::cerr << ind << "Class_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <<" << selector->ns_name() << ">>";
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << std::endl;
  } else if (Cast<Id_Selector>(node)) {
    Id_Selector_Ptr selector = Cast<Id_Selector>(node);
    std::cerr << ind << "Id_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <<" << selector->ns_name() << ">>";
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << std::endl;
  } else if (Cast<Type_Selector>(node)) {
    Type_Selector_Ptr selector = Cast<Type_Selector>(node);
    std::cerr << ind << "Type_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <" << selector->hash() << ">";
    std::cerr << " <<" << selector->ns_name() << ">>";
    std::cerr << (selector->is_optional() ? " [is_optional]": " -");
    std::cerr << (selector->has_parent_ref() ? " [has-parent]": " -");
    std::cerr << (selector->has_line_break() ? " [line-break]": " -");
    std::cerr << (selector->has_line_feed() ? " [line-feed]": " -");
    std::cerr << " <" << prettyprint(selector->pstate().token.ws_before()) << ">";
    std::cerr << std::endl;
  } else if (Cast<Placeholder_Selector>(node)) {

    Placeholder_Selector_Ptr selector = Cast<Placeholder_Selector>(node);
    std::cerr << ind << "Placeholder_Selector [" << selector->ns_name() << "] " << selector;
    std::cerr << " (" << pstate_source_position(selector) << ")"
      << " <" << selector->hash() << ">"
      << " [@media:" << selector->media_block() << "]"
      << (selector->is_optional() ? " [is_optional]": " -")
      << (selector->has_line_break() ? " [line-break]": " -")
      << (selector->has_line_feed() ? " [line-feed]": " -")
    << std::endl;

  } else if (Cast<Simple_Selector>(node)) {
    Simple_Selector* selector = Cast<Simple_Selector>(node);
    std::cerr << ind << "Simple_Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << (selector->has_line_break() ? " [line-break]": " -") << (selector->has_line_feed() ? " [line-feed]": " -") << std::endl;

  } else if (Cast<Selector_Schema>(node)) {
    Selector_Schema_Ptr selector = Cast<Selector_Schema>(node);
    std::cerr << ind << "Selector_Schema " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")"
      << " [@media:" << selector->media_block() << "]"
      << (selector->connect_parent() ? " [connect-parent]": " -")
    << std::endl;

    debug_ast(selector->contents(), ind + " ");
    // for(auto i : selector->elements()) { debug_ast(i, ind + " ", env); }

  } else if (Cast<Selector>(node)) {
    Selector_Ptr selector = Cast<Selector>(node);
    std::cerr << ind << "Selector " << selector;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << (selector->has_line_break() ? " [line-break]": " -")
      << (selector->has_line_feed() ? " [line-feed]": " -")
    << std::endl;

  } else if (Cast<Media_Query_Expression>(node)) {
    Media_Query_Expression_Ptr block = Cast<Media_Query_Expression>(node);
    std::cerr << ind << "Media_Query_Expression " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << (block->is_interpolated() ? " [is_interpolated]": " -")
    << std::endl;
    debug_ast(block->feature(), ind + " feature) ");
    debug_ast(block->value(), ind + " value) ");

  } else if (Cast<Media_Query>(node)) {
    Media_Query_Ptr block = Cast<Media_Query>(node);
    std::cerr << ind << "Media_Query " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << (block->is_negated() ? " [is_negated]": " -")
      << (block->is_restricted() ? " [is_restricted]": " -")
    << std::endl;
    debug_ast(block->media_type(), ind + " ");
    for(const auto& i : block->elements()) { debug_ast(i, ind + " ", env); }

  } else if (Cast<Media_Block>(node)) {
    Media_Block_Ptr block = Cast<Media_Block>(node);
    std::cerr << ind << "Media_Block " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->media_queries(), ind + " =@ ");
    if (block->block()) for(const Statement_Obj& i : block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Supports_Block>(node)) {
    Supports_Block_Ptr block = Cast<Supports_Block>(node);
    std::cerr << ind << "Supports_Block " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->condition(), ind + " =@ ");
    debug_ast(block->block(), ind + " <>");
  } else if (Cast<Supports_Operator>(node)) {
    Supports_Operator_Ptr block = Cast<Supports_Operator>(node);
    std::cerr << ind << "Supports_Operator " << block;
    std::cerr << " (" << pstate_source_position(node) << ")"
    << std::endl;
    debug_ast(block->left(), ind + " left) ");
    debug_ast(block->right(), ind + " right) ");
  } else if (Cast<Supports_Negation>(node)) {
    Supports_Negation_Ptr block = Cast<Supports_Negation>(node);
    std::cerr << ind << "Supports_Negation " << block;
    std::cerr << " (" << pstate_source_position(node) << ")"
    << std::endl;
    debug_ast(block->condition(), ind + " condition) ");
  } else if (Cast<At_Root_Query>(node)) {
    At_Root_Query_Ptr block = Cast<At_Root_Query>(node);
    std::cerr << ind << "At_Root_Query " << block;
    std::cerr << " (" << pstate_source_position(node) << ")"
    << std::endl;
    debug_ast(block->feature(), ind + " feature) ");
    debug_ast(block->value(), ind + " value) ");
  } else if (Cast<Supports_Declaration>(node)) {
    Supports_Declaration_Ptr block = Cast<Supports_Declaration>(node);
    std::cerr << ind << "Supports_Declaration " << block;
    std::cerr << " (" << pstate_source_position(node) << ")"
    << std::endl;
    debug_ast(block->feature(), ind + " feature) ");
    debug_ast(block->value(), ind + " value) ");
  } else if (Cast<Block>(node)) {
    Block_Ptr root_block = Cast<Block>(node);
    std::cerr << ind << "Block " << root_block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    if (root_block->is_root()) std::cerr << " [root]";
    std::cerr << " " << root_block->tabs() << std::endl;
    for(const Statement_Obj& i : root_block->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Warning>(node)) {
    Warning_Ptr block = Cast<Warning>(node);
    std::cerr << ind << "Warning " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->message(), ind + " : ");
  } else if (Cast<Error>(node)) {
    Error_Ptr block = Cast<Error>(node);
    std::cerr << ind << "Error " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
  } else if (Cast<Debug>(node)) {
    Debug_Ptr block = Cast<Debug>(node);
    std::cerr << ind << "Debug " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->value(), ind + " ");
  } else if (Cast<Comment>(node)) {
    Comment_Ptr block = Cast<Comment>(node);
    std::cerr << ind << "Comment " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() <<
      " <" << prettyprint(block->pstate().token.ws_before()) << ">" << std::endl;
    debug_ast(block->text(), ind + "// ", env);
  } else if (Cast<If>(node)) {
    If_Ptr block = Cast<If>(node);
    std::cerr << ind << "If " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->predicate(), ind + " = ");
    debug_ast(block->block(), ind + " <>");
    debug_ast(block->alternative(), ind + " ><");
  } else if (Cast<Return>(node)) {
    Return_Ptr block = Cast<Return>(node);
    std::cerr << ind << "Return " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
  } else if (Cast<Extension>(node)) {
    Extension_Ptr block = Cast<Extension>(node);
    std::cerr << ind << "Extension " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->selector(), ind + "-> ", env);
  } else if (Cast<Content>(node)) {
    Content_Ptr block = Cast<Content>(node);
    std::cerr << ind << "Content " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->arguments(), ind + " args: ", env);
  } else if (Cast<Import_Stub>(node)) {
    Import_Stub_Ptr block = Cast<Import_Stub>(node);
    std::cerr << ind << "Import_Stub " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << block->imp_path() << "] ";
    std::cerr << " " << block->tabs() << std::endl;
  } else if (Cast<Import>(node)) {
    Import_Ptr block = Cast<Import>(node);
    std::cerr << ind << "Import " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    // std::vector<std::string>         files_;
    for (auto imp : block->urls()) debug_ast(imp, ind + "@: ", env);
    debug_ast(block->import_queries(), ind + "@@ ");
  } else if (Cast<Assignment>(node)) {
    Assignment_Ptr block = Cast<Assignment>(node);
    std::cerr << ind << "Assignment " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " <<" << block->variable() << ">> " << block->tabs() << std::endl;
    debug_ast(block->value(), ind + "=", env);
  } else if (Cast<Declaration>(node)) {
    Declaration_Ptr block = Cast<Declaration>(node);
    std::cerr << ind << "Declaration " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [is_custom_property: " << block->is_custom_property() << "] ";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->property(), ind + " prop: ", env);
    debug_ast(block->value(), ind + " value: ", env);
    debug_ast(block->block(), ind + " ", env);
  } else if (Cast<Keyframe_Rule>(node)) {
    Keyframe_Rule_Ptr has_block = Cast<Keyframe_Rule>(node);
    std::cerr << ind << "Keyframe_Rule " << has_block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << has_block->tabs() << std::endl;
    if (has_block->name()) debug_ast(has_block->name(), ind + "@");
    if (has_block->block()) for(const Statement_Obj& i : has_block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Directive>(node)) {
    Directive_Ptr block = Cast<Directive>(node);
    std::cerr << ind << "Directive " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << block->keyword() << "] " << block->tabs() << std::endl;
    debug_ast(block->selector(), ind + "~", env);
    debug_ast(block->value(), ind + "+", env);
    if (block->block()) for(const Statement_Obj& i : block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Each>(node)) {
    Each_Ptr block = Cast<Each>(node);
    std::cerr << ind << "Each " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    if (block->block()) for(const Statement_Obj& i : block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<For>(node)) {
    For_Ptr block = Cast<For>(node);
    std::cerr << ind << "For " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    if (block->block()) for(const Statement_Obj& i : block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<While>(node)) {
    While_Ptr block = Cast<While>(node);
    std::cerr << ind << "While " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << block->tabs() << std::endl;
    if (block->block()) for(const Statement_Obj& i : block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Definition>(node)) {
    Definition_Ptr block = Cast<Definition>(node);
    std::cerr << ind << "Definition " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [name: " << block->name() << "] ";
    std::cerr << " [type: " << (block->type() == Sass::Definition::Type::MIXIN ? "Mixin " : "Function ") << "] ";
    // this seems to lead to segfaults some times?
    // std::cerr << " [signature: " << block->signature() << "] ";
    std::cerr << " [native: " << block->native_function() << "] ";
    std::cerr << " " << block->tabs() << std::endl;
    debug_ast(block->parameters(), ind + " params: ", env);
    if (block->block()) debug_ast(block->block(), ind + " ", env);
  } else if (Cast<Mixin_Call>(node)) {
    Mixin_Call_Ptr block = Cast<Mixin_Call>(node);
    std::cerr << ind << "Mixin_Call " << block << " " << block->tabs();
    std::cerr << " (" << pstate_source_position(block) << ")";
    std::cerr << " [" <<  block->name() << "]";
    std::cerr << " [has_content: " << block->has_content() << "] " << std::endl;
    debug_ast(block->arguments(), ind + " args: ", env);
    debug_ast(block->block_parameters(), ind + " block_params: ", env);
    if (block->block()) debug_ast(block->block(), ind + " ", env);
  } else if (Ruleset_Ptr ruleset = Cast<Ruleset>(node)) {
    std::cerr << ind << "Ruleset " << ruleset;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [indent: " << ruleset->tabs() << "]";
    std::cerr << (ruleset->is_invisible() ? " [INVISIBLE]" : "");
    std::cerr << (ruleset->is_root() ? " [root]" : "");
    std::cerr << std::endl;
    debug_ast(ruleset->selector(), ind + ">");
    debug_ast(ruleset->block(), ind + " ");
  } else if (Cast<Block>(node)) {
    Block_Ptr block = Cast<Block>(node);
    std::cerr << ind << "Block " << block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << (block->is_invisible() ? " [INVISIBLE]" : "");
    std::cerr << " [indent: " << block->tabs() << "]" << std::endl;
    for(const Statement_Obj& i : block->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Variable>(node)) {
    Variable_Ptr expression = Cast<Variable>(node);
    std::cerr << ind << "Variable " << expression;
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << expression->name() << "]" << std::endl;
    std::string name(expression->name());
    if (env && env->has(name)) debug_ast(Cast<Expression>((*env)[name]), ind + " -> ", env);
  } else if (Cast<Function_Call>(node)) {
    Function_Call_Ptr expression = Cast<Function_Call>(node);
    std::cerr << ind << "Function_Call " << expression;
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << expression->name() << "]";
    if (expression->is_delayed()) std::cerr << " [delayed]";
    if (expression->is_interpolant()) std::cerr << " [interpolant]";
    if (expression->is_css()) std::cerr << " [css]";
    std::cerr << std::endl;
    debug_ast(expression->arguments(), ind + " args: ", env);
    debug_ast(expression->func(), ind + " func: ", env);
  } else if (Cast<Function>(node)) {
    Function_Ptr expression = Cast<Function>(node);
    std::cerr << ind << "Function " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    if (expression->is_css()) std::cerr << " [css]";
    std::cerr << std::endl;
    debug_ast(expression->definition(), ind + " definition: ", env);
  } else if (Cast<Arguments>(node)) {
    Arguments_Ptr expression = Cast<Arguments>(node);
    std::cerr << ind << "Arguments " << expression;
    if (expression->is_delayed()) std::cerr << " [delayed]";
    std::cerr << " (" << pstate_source_position(node) << ")";
    if (expression->has_named_arguments()) std::cerr << " [has_named_arguments]";
    if (expression->has_rest_argument()) std::cerr << " [has_rest_argument]";
    if (expression->has_keyword_argument()) std::cerr << " [has_keyword_argument]";
    std::cerr << std::endl;
    for(const Argument_Obj& i : expression->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Argument>(node)) {
    Argument_Ptr expression = Cast<Argument>(node);
    std::cerr << ind << "Argument " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << expression->value().ptr() << "]";
    std::cerr << " [name: " << expression->name() << "] ";
    std::cerr << " [rest: " << expression->is_rest_argument() << "] ";
    std::cerr << " [keyword: " << expression->is_keyword_argument() << "] " << std::endl;
    debug_ast(expression->value(), ind + " value: ", env);
  } else if (Cast<Parameters>(node)) {
    Parameters_Ptr expression = Cast<Parameters>(node);
    std::cerr << ind << "Parameters " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [has_optional: " << expression->has_optional_parameters() << "] ";
    std::cerr << " [has_rest: " << expression->has_rest_parameter() << "] ";
    std::cerr << std::endl;
    for(const Parameter_Obj& i : expression->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Parameter>(node)) {
    Parameter_Ptr expression = Cast<Parameter>(node);
    std::cerr << ind << "Parameter " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [name: " << expression->name() << "] ";
    std::cerr << " [default: " << expression->default_value().ptr() << "] ";
    std::cerr << " [rest: " << expression->is_rest_parameter() << "] " << std::endl;
  } else if (Cast<Unary_Expression>(node)) {
    Unary_Expression_Ptr expression = Cast<Unary_Expression>(node);
    std::cerr << ind << "Unary_Expression " << expression;
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " [delayed: " << expression->is_delayed() << "] ";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << expression->type() << "]" << std::endl;
    debug_ast(expression->operand(), ind + " operand: ", env);
  } else if (Cast<Binary_Expression>(node)) {
    Binary_Expression_Ptr expression = Cast<Binary_Expression>(node);
    std::cerr << ind << "Binary_Expression " << expression;
    if (expression->is_interpolant()) std::cerr << " [is interpolant] ";
    if (expression->is_left_interpolant()) std::cerr << " [left interpolant] ";
    if (expression->is_right_interpolant()) std::cerr << " [right interpolant] ";
    std::cerr << " [delayed: " << expression->is_delayed() << "] ";
    std::cerr << " [ws_before: " << expression->op().ws_before << "] ";
    std::cerr << " [ws_after: " << expression->op().ws_after << "] ";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << expression->type_name() << "]" << std::endl;
    debug_ast(expression->left(), ind + " left:  ", env);
    debug_ast(expression->right(), ind + " right: ", env);
  } else if (Cast<Map>(node)) {
    Map_Ptr expression = Cast<Map>(node);
    std::cerr << ind << "Map " << expression;
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [Hashed]" << std::endl;
    for (const auto& i : expression->elements()) {
      debug_ast(i.first, ind + " key: ");
      debug_ast(i.second, ind + " val: ");
    }
  } else if (Cast<List>(node)) {
    List_Ptr expression = Cast<List>(node);
    std::cerr << ind << "List " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " (" << expression->length() << ") " <<
      (expression->separator() == SASS_COMMA ? "Comma " : expression->separator() == SASS_HASH ? "Map " : "Space ") <<
      " [delayed: " << expression->is_delayed() << "] " <<
      " [interpolant: " << expression->is_interpolant() << "] " <<
      " [listized: " << expression->from_selector() << "] " <<
      " [arglist: " << expression->is_arglist() << "] " <<
      " [bracketed: " << expression->is_bracketed() << "] " <<
      " [expanded: " << expression->is_expanded() << "] " <<
      " [hash: " << expression->hash() << "] " <<
      std::endl;
    for(const auto& i : expression->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Boolean>(node)) {
    Boolean_Ptr expression = Cast<Boolean>(node);
    std::cerr << ind << "Boolean " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " [" << expression->value() << "]" << std::endl;
  } else if (Cast<Color_RGBA>(node)) {
    Color_RGBA_Ptr expression = Cast<Color_RGBA>(node);
    std::cerr << ind << "Color " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [name: " << expression->disp() << "] ";
    std::cerr << " [delayed: " << expression->is_delayed() << "] ";
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " rgba[" << expression->r() << ":"  << expression->g() << ":" << expression->b() << "@" << expression->a() << "]" << std::endl;
  } else if (Cast<Color_HSLA>(node)) {
    Color_HSLA_Ptr expression = Cast<Color_HSLA>(node);
    std::cerr << ind << "Color " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [name: " << expression->disp() << "] ";
    std::cerr << " [delayed: " << expression->is_delayed() << "] ";
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " hsla[" << expression->h() << ":"  << expression->s() << ":" << expression->l() << "@" << expression->a() << "]" << std::endl;
  } else if (Cast<Number>(node)) {
    Number_Ptr expression = Cast<Number>(node);
    std::cerr << ind << "Number " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [delayed: " << expression->is_delayed() << "] ";
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] ";
    std::cerr << " [" << expression->value() << expression->unit() << "]" <<
      " [hash: " << expression->hash() << "] " <<
      std::endl;
  } else if (Cast<Null>(node)) {
    Null_Ptr expression = Cast<Null>(node);
    std::cerr << ind << "Null " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [interpolant: " << expression->is_interpolant() << "] "
      // " [hash: " << expression->hash() << "] "
      << std::endl;
  } else if (Cast<String_Quoted>(node)) {
    String_Quoted_Ptr expression = Cast<String_Quoted>(node);
    std::cerr << ind << "String_Quoted " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << prettyprint(expression->value()) << "]";
    if (expression->is_delayed()) std::cerr << " [delayed]";
    if (expression->is_interpolant()) std::cerr << " [interpolant]";
    if (expression->quote_mark()) std::cerr << " [quote_mark: " << expression->quote_mark() << "]";
    std::cerr << " <" << prettyprint(expression->pstate().token.ws_before()) << ">" << std::endl;
  } else if (Cast<String_Constant>(node)) {
    String_Constant_Ptr expression = Cast<String_Constant>(node);
    std::cerr << ind << "String_Constant " << expression;
    if (expression->concrete_type()) {
      std::cerr << " " << expression->concrete_type();
    }
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " [" << prettyprint(expression->value()) << "]";
    if (expression->is_delayed()) std::cerr << " [delayed]";
    if (expression->is_interpolant()) std::cerr << " [interpolant]";
    std::cerr << " <" << prettyprint(expression->pstate().token.ws_before()) << ">" << std::endl;
  } else if (Cast<String_Schema>(node)) {
    String_Schema_Ptr expression = Cast<String_Schema>(node);
    std::cerr << ind << "String_Schema " << expression;
    std::cerr << " (" << pstate_source_position(expression) << ")";
    std::cerr << " " << expression->concrete_type();
    std::cerr << " (" << pstate_source_position(node) << ")";
    if (expression->css()) std::cerr << " [css]";
    if (expression->is_delayed()) std::cerr << " [delayed]";
    if (expression->is_interpolant()) std::cerr << " [is interpolant]";
    if (expression->has_interpolant()) std::cerr << " [has interpolant]";
    if (expression->is_left_interpolant()) std::cerr << " [left interpolant] ";
    if (expression->is_right_interpolant()) std::cerr << " [right interpolant] ";
    std::cerr << " <" << prettyprint(expression->pstate().token.ws_before()) << ">" << std::endl;
    for(const auto& i : expression->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<String>(node)) {
    String_Ptr expression = Cast<String>(node);
    std::cerr << ind << "String " << expression;
    std::cerr << " " << expression->concrete_type();
    std::cerr << " (" << pstate_source_position(node) << ")";
    if (expression->is_interpolant()) std::cerr << " [interpolant]";
    std::cerr << " <" << prettyprint(expression->pstate().token.ws_before()) << ">" << std::endl;
  } else if (Cast<Expression>(node)) {
    Expression_Ptr expression = Cast<Expression>(node);
    std::cerr << ind << "Expression " << expression;
    std::cerr << " (" << pstate_source_position(node) << ")";
    switch (expression->concrete_type()) {
      case Expression::Type::NONE: std::cerr << " [NONE]"; break;
      case Expression::Type::BOOLEAN: std::cerr << " [BOOLEAN]"; break;
      case Expression::Type::NUMBER: std::cerr << " [NUMBER]"; break;
      case Expression::Type::COLOR: std::cerr << " [COLOR]"; break;
      case Expression::Type::STRING: std::cerr << " [STRING]"; break;
      case Expression::Type::LIST: std::cerr << " [LIST]"; break;
      case Expression::Type::MAP: std::cerr << " [MAP]"; break;
      case Expression::Type::SELECTOR: std::cerr << " [SELECTOR]"; break;
      case Expression::Type::NULL_VAL: std::cerr << " [NULL_VAL]"; break;
      case Expression::Type::C_WARNING: std::cerr << " [C_WARNING]"; break;
      case Expression::Type::C_ERROR: std::cerr << " [C_ERROR]"; break;
      case Expression::Type::FUNCTION: std::cerr << " [FUNCTION]"; break;
      case Expression::Type::NUM_TYPES: std::cerr << " [NUM_TYPES]"; break;
      case Expression::Type::VARIABLE: std::cerr << " [VARIABLE]"; break;
      case Expression::Type::FUNCTION_VAL: std::cerr << " [FUNCTION_VAL]"; break;
      case Expression::Type::PARENT: std::cerr << " [PARENT]"; break;
    }
    std::cerr << std::endl;
  } else if (Cast<Has_Block>(node)) {
    Has_Block_Ptr has_block = Cast<Has_Block>(node);
    std::cerr << ind << "Has_Block " << has_block;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << has_block->tabs() << std::endl;
    if (has_block->block()) for(const Statement_Obj& i : has_block->block()->elements()) { debug_ast(i, ind + " ", env); }
  } else if (Cast<Statement>(node)) {
    Statement_Ptr statement = Cast<Statement>(node);
    std::cerr << ind << "Statement " << statement;
    std::cerr << " (" << pstate_source_position(node) << ")";
    std::cerr << " " << statement->tabs() << std::endl;
  }

  if (ind == "") std::cerr << "####################################################################\n";
}

inline void debug_node(Node* node, std::string ind = "")
{
  if (ind == "") std::cerr << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  if (node->isCombinator()) {
    std::cerr << ind;
    std::cerr << "Combinator ";
    std::cerr << node << " ";
    if (node->got_line_feed) std::cerr << "[LF] ";
    switch (node->combinator()) {
      case Complex_Selector::ADJACENT_TO: std::cerr << "{+} "; break;
      case Complex_Selector::PARENT_OF:   std::cerr << "{>} "; break;
      case Complex_Selector::PRECEDES:    std::cerr << "{~} "; break;
      case Complex_Selector::REFERENCE:   std::cerr << "{@} "; break;
      case Complex_Selector::ANCESTOR_OF: std::cerr << "{ } "; break;
    }
    std::cerr << std::endl;
    // debug_ast(node->combinator(), ind + "  ");
  } else if (node->isSelector()) {
    std::cerr << ind;
    std::cerr << "Selector ";
    std::cerr << node << " ";
    if (node->got_line_feed) std::cerr << "[LF] ";
    std::cerr << std::endl;
    debug_ast(node->selector(), ind + "  ");
  } else if (node->isCollection()) {
    std::cerr << ind;
    std::cerr << "Collection ";
    std::cerr << node << " ";
    if (node->got_line_feed) std::cerr << "[LF] ";
    std::cerr << std::endl;
    for(auto n : (*node->collection())) {
      debug_node(&n, ind + "  ");
    }
  } else if (node->isNil()) {
    std::cerr << ind;
    std::cerr << "Nil ";
    std::cerr << node << " ";
    if (node->got_line_feed) std::cerr << "[LF] ";
    std::cerr << std::endl;
  } else {
    std::cerr << ind;
    std::cerr << "OTHER ";
    std::cerr << node << " ";
    if (node->got_line_feed) std::cerr << "[LF] ";
    std::cerr << std::endl;
  }
  if (ind == "") std::cerr << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
}

/*
inline void debug_ast(const AST_Node_Ptr node, std::string ind = "", Env* env = 0)
{
  debug_ast(const_cast<AST_Node_Ptr>(node), ind, env);
}
*/
inline void debug_node(const Node* node, std::string ind = "")
{
  debug_node(const_cast<Node*>(node), ind);
}

inline void debug_subset_map(Sass::Subset_Map& map, std::string ind = "")
{
  if (ind == "") std::cerr << "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  for(auto const &it : map.values()) {
    debug_ast(it.first, ind + "first: ");
    debug_ast(it.second, ind + "second: ");
  }
  if (ind == "") std::cerr << "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
}

inline void debug_subset_entries(SubSetMapPairs* entries, std::string ind = "")
{
  if (ind == "") std::cerr << "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  for(auto const &pair : *entries) {
    debug_ast(pair.first, ind + "first: ");
    debug_ast(pair.second, ind + "second: ");
  }
  if (ind == "") std::cerr << "#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
}

#endif // SASS_DEBUGGER
#ifndef SASS_AST2C_H
#define SASS_AST2C_H


namespace Sass {

  class AST2C : public Operation_CRTP<union Sass_Value*, AST2C> {

  public:

    AST2C() { }
    ~AST2C() { }

    union Sass_Value* operator()(Boolean_Ptr);
    union Sass_Value* operator()(Number_Ptr);
    union Sass_Value* operator()(Color_RGBA_Ptr);
    union Sass_Value* operator()(Color_HSLA_Ptr);
    union Sass_Value* operator()(String_Constant_Ptr);
    union Sass_Value* operator()(String_Quoted_Ptr);
    union Sass_Value* operator()(Custom_Warning_Ptr);
    union Sass_Value* operator()(Custom_Error_Ptr);
    union Sass_Value* operator()(List_Ptr);
    union Sass_Value* operator()(Map_Ptr);
    union Sass_Value* operator()(Null_Ptr);
    union Sass_Value* operator()(Arguments_Ptr);
    union Sass_Value* operator()(Argument_Ptr);

    // return sass error if type is not supported
    union Sass_Value* fallback(AST_Node_Ptr x)
    { return sass_make_error("unknown type for C-API"); }

  };

}

#endif
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

#ifndef CCAN_JSON_H
#define CCAN_JSON_H


typedef enum {
  JSON_NULL,
  JSON_BOOL,
  JSON_STRING,
  JSON_NUMBER,
  JSON_ARRAY,
  JSON_OBJECT,
} JsonTag;

typedef struct JsonNode JsonNode;

struct JsonNode
{
  /* only if parent is an object or array (NULL otherwise) */
  JsonNode *parent;
  JsonNode *prev, *next;

  /* only if parent is an object (NULL otherwise) */
  char *key; /* Must be valid UTF-8. */

  JsonTag tag;
  union {
    /* JSON_BOOL */
    bool bool_;

    /* JSON_STRING */
    char *string_; /* Must be valid UTF-8. */

    /* JSON_NUMBER */
    double number_;

    /* JSON_ARRAY */
    /* JSON_OBJECT */
    struct {
      JsonNode *head, *tail;
    } children;
  };
};

/*** Encoding, decoding, and validation ***/

JsonNode   *json_decode         (const char *json);
char       *json_encode         (const JsonNode *node);
char       *json_encode_string  (const char *str);
char       *json_stringify      (const JsonNode *node, const char *space);
void        json_delete         (JsonNode *node);

bool        json_validate       (const char *json);

/*** Lookup and traversal ***/

JsonNode   *json_find_element   (JsonNode *array, int index);
JsonNode   *json_find_member    (JsonNode *object, const char *key);

JsonNode   *json_first_child    (const JsonNode *node);

#define json_foreach(i, object_or_array)            \
  for ((i) = json_first_child(object_or_array);   \
     (i) != NULL;                               \
     (i) = (i)->next)

/*** Construction and manipulation ***/

JsonNode *json_mknull(void);
JsonNode *json_mkbool(bool b);
JsonNode *json_mkstring(const char *s);
JsonNode *json_mknumber(double n);
JsonNode *json_mkarray(void);
JsonNode *json_mkobject(void);

void json_append_element(JsonNode *array, JsonNode *element);
void json_prepend_element(JsonNode *array, JsonNode *element);
void json_append_member(JsonNode *object, const char *key, JsonNode *value);
void json_prepend_member(JsonNode *object, const char *key, JsonNode *value);

void json_remove_from_parent(JsonNode *node);

/*** Debugging ***/

/*
 * Look for structure and encoding problems in a JsonNode or its descendents.
 *
 * If a problem is detected, return false, writing a description of the problem
 * to errmsg (unless errmsg is NULL).
 */
bool json_check(const JsonNode *node, char errmsg[256]);

#endif
#ifndef SASS_REMOVE_PLACEHOLDERS_H
#define SASS_REMOVE_PLACEHOLDERS_H



namespace Sass {


    class Remove_Placeholders : public Operation_CRTP<void, Remove_Placeholders> {

    public:
      Selector_List_Ptr remove_placeholders(Selector_List_Ptr);

    public:
        Remove_Placeholders();
        ~Remove_Placeholders() { }

        void operator()(Block_Ptr);
        void operator()(Ruleset_Ptr);
        void operator()(Media_Block_Ptr);
        void operator()(Supports_Block_Ptr);
        void operator()(Directive_Ptr);

      // ignore missed types
        template <typename U>
      void fallback(U x) {}
    };

}

#endif
#ifndef SASS_FN_MISCS_H
#define SASS_FN_MISCS_H


namespace Sass {

  namespace Functions {

    extern Signature type_of_sig;
    extern Signature variable_exists_sig;
    extern Signature global_variable_exists_sig;
    extern Signature function_exists_sig;
    extern Signature mixin_exists_sig;
    extern Signature feature_exists_sig;
    extern Signature call_sig;
    extern Signature not_sig;
    extern Signature if_sig;
    extern Signature set_nth_sig;
    extern Signature content_exists_sig;
    extern Signature get_function_sig;

    BUILT_IN(type_of);
    BUILT_IN(variable_exists);
    BUILT_IN(global_variable_exists);
    BUILT_IN(function_exists);
    BUILT_IN(mixin_exists);
    BUILT_IN(feature_exists);
    BUILT_IN(call);
    BUILT_IN(sass_not);
    BUILT_IN(sass_if);
    BUILT_IN(set_nth);
    BUILT_IN(content_exists);
    BUILT_IN(get_function);

  }

}

#endif
#ifndef SASS_FN_NUMBERS_H
#define SASS_FN_NUMBERS_H


namespace Sass {

  namespace Functions {

    // return a number object (copied since we want to have reduced units)
    #define ARGN(argname) get_arg_n(argname, env, sig, pstate, traces) // Number copy

    extern Signature percentage_sig;
    extern Signature round_sig;
    extern Signature ceil_sig;
    extern Signature floor_sig;
    extern Signature abs_sig;
    extern Signature min_sig;
    extern Signature max_sig;
    extern Signature inspect_sig;
    extern Signature random_sig;
    extern Signature unique_id_sig;
    extern Signature unit_sig;
    extern Signature unitless_sig;
    extern Signature comparable_sig;

    BUILT_IN(percentage);
    BUILT_IN(round);
    BUILT_IN(ceil);
    BUILT_IN(floor);
    BUILT_IN(abs);
    BUILT_IN(min);
    BUILT_IN(max);
    BUILT_IN(inspect);
    BUILT_IN(random);
    BUILT_IN(unique_id);
    BUILT_IN(unit);
    BUILT_IN(unitless);
    BUILT_IN(comparable);

  }

}

#endif
#ifndef SASS_SASS_VALUES_H
#define SASS_SASS_VALUES_H


struct Sass_Unknown {
  enum Sass_Tag tag;
};

struct Sass_Boolean {
  enum Sass_Tag tag;
  bool          value;
};

struct Sass_Number {
  enum Sass_Tag tag;
  double        value;
  char*         unit;
};

struct Sass_Color {
  enum Sass_Tag tag;
  double        r;
  double        g;
  double        b;
  double        a;
};

struct Sass_String {
  enum Sass_Tag tag;
  bool          quoted;
  char*         value;
};

struct Sass_List {
  enum Sass_Tag       tag;
  enum Sass_Separator separator;
  bool                is_bracketed;
  size_t              length;
  // null terminated "array"
  union Sass_Value**  values;
};

struct Sass_Map {
  enum Sass_Tag        tag;
  size_t               length;
  struct Sass_MapPair* pairs;
};

struct Sass_Null {
  enum Sass_Tag tag;
};

struct Sass_Error {
  enum Sass_Tag tag;
  char*         message;
};

struct Sass_Warning {
  enum Sass_Tag tag;
  char*         message;
};

union Sass_Value {
  struct Sass_Unknown unknown;
  struct Sass_Boolean boolean;
  struct Sass_Number  number;
  struct Sass_Color   color;
  struct Sass_String  string;
  struct Sass_List    list;
  struct Sass_Map     map;
  struct Sass_Null    null;
  struct Sass_Error   error;
  struct Sass_Warning warning;
};

struct Sass_MapPair {
  union Sass_Value* key;
  union Sass_Value* value;
};

#endif
#ifndef SASS_FN_STRINGS_H
#define SASS_FN_STRINGS_H


namespace Sass {

  namespace Functions {

    extern Signature unquote_sig;
    extern Signature quote_sig;
    extern Signature str_length_sig;
    extern Signature str_insert_sig;
    extern Signature str_index_sig;
    extern Signature str_slice_sig;
    extern Signature to_upper_case_sig;
    extern Signature to_lower_case_sig;
    extern Signature length_sig;

    BUILT_IN(sass_unquote);
    BUILT_IN(sass_quote);
    BUILT_IN(str_length);
    BUILT_IN(str_insert);
    BUILT_IN(str_index);
    BUILT_IN(str_slice);
    BUILT_IN(to_upper_case);
    BUILT_IN(to_lower_case);
    BUILT_IN(length);

  }

}

#endif
#ifndef SASS_BIND_H
#define SASS_BIND_H


namespace Sass {

  void bind(std::string type, std::string name, Parameters_Obj, Arguments_Obj, Env*, Eval*, Backtraces& traces);

}

#endif
#ifndef SASS_FN_MAPS_H
#define SASS_FN_MAPS_H


namespace Sass {

  namespace Functions {

    #define ARGM(argname, argtype) get_arg_m(argname, env, sig, pstate, traces)

    extern Signature map_get_sig;
    extern Signature map_merge_sig;
    extern Signature map_remove_sig;
    extern Signature map_keys_sig;
    extern Signature map_values_sig;
    extern Signature map_has_key_sig;

    BUILT_IN(map_get);
    BUILT_IN(map_merge);
    BUILT_IN(map_remove);
    BUILT_IN(map_keys);
    BUILT_IN(map_values);
    BUILT_IN(map_has_key);

  }

}

#endif
