//! \file ConfBlock.cpp
//! Brief description

//! \date   Jan 17, 2009
//! \author fuan
//!
//! Detailed description

#include <iostream>
#include <fstream>
#include <boost/algorithm/string_regex.hpp>

#include "ConfBlock.h"

using boost::lexical_cast;

namespace confparser {

// ============================ //
// Constructors and destructors //
// ============================ //

ConfBlock::ConfBlock( std::string blockName, int level, ConfBlock* parent )
    : blockName_( blockName ),
      level_( level ),
      props_(),
      child_( NULL ),
      sibling_( NULL ),
      parent_( parent ) {}

ConfBlock::~ConfBlock() {

  if ( child_ != NULL ) delete child_;
  if ( sibling_ != NULL ) delete sibling_;

}

void ConfBlock::unlink() {

  // Remove the block from the siblings unless it is the outermost block
  if ( parent_ != NULL ) {

    // If the block is the first child, then its first sibling becomes the first
    // child of the parent
    if ( parent_->child_ == this ) {
      parent_->child_ = sibling_;
    // Else find the previous sibling and unlink the block from the chain of
    // siblings
    } else {
      ConfBlock* prevSib = parent_->child_;
      while ( prevSib->sibling_ != this ) prevSib = prevSib->sibling_;
      prevSib->sibling_ = sibling_;
    }
  }

  // Then set the sibling to NULL in order to save it from destruction and call
  // destructor
  sibling_ = NULL;
  delete this;
}

// ======= //
// Getters //
// ======= //

std::string ConfBlock::getQualifiedName() {

  std::string qName = blockName_;
  ConfBlock* searchBlock = this;

  while ( searchBlock->level_ > 0 ) {
    searchBlock = searchBlock->parent_;
    qName = searchBlock->blockName_ + "." + qName;
  }

  return qName;
}

ConfBlock* ConfBlock::findSibling( std::string name ) {

  ConfBlock* sib = this;

  while ( sib->sibling_ != NULL ) {
    sib = sib->sibling_;
    if ( sib->blockName_ == name ) return sib;
  }

  return NULL;
}

ConfBlock* ConfBlock::findRec( std::string name ) {

  ConfBlock* sib = this;
  int initLvl = level_;

  while ( sib->blockName_ != name ) {
    if ( sib->child_ != NULL ) sib = sib->child_;
    else if ( sib->sibling_ != NULL ) sib = sib->sibling_;
    else if ( sib->level_ > initLvl && sib->level_ > 0 ) sib = sib->parent_;
    else return NULL;
  }

  return sib;
}

// ======= //
// Setters //
// ======= //

void ConfBlock::writeConfigFile( std::string fileName ) {

  std::ofstream f( fileName.c_str() );
  f << "# Automatically generated config file\n\n";
  writeConfigFileRec( f, level_ );
  f.close();

}

// ========================== //
// Protected helper functions //
// ========================== //

void ConfBlock::writeConfigFileRec( std::ofstream &fileHandle, int initLvl ) {

  std::string indent( level_ == initLvl ? 0 : 2 * ( level_ - initLvl - 1), ' ' );

  // Start new block unless we are the outermost
  if ( level_ > initLvl ) fileHandle << indent << blockName_ << " {\n\n";

  // Write out all data of current block
  std::map<std::string,std::string>::iterator mit;
  for ( mit = props_.begin(); mit != props_.end(); ++mit ) {
    fileHandle << indent << ( level_ > initLvl ? "  " : "" ) << mit->first << " " << mit->second << ";\n";
  }
  fileHandle << '\n';

  // Process children first, if any
  if ( child_ != NULL ) child_->writeConfigFileRec( fileHandle, initLvl );

  // End block unless we are the outermost
  if ( level_ > initLvl ) fileHandle << indent << "}\n\n";

  // Continue with next sibling if any
  if ( sibling_ != NULL ) sibling_->writeConfigFileRec( fileHandle, initLvl );

}

} // namespace confparser
