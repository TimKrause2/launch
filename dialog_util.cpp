#include "dialog_util.h"
#include <ios>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <cmath>


#define RAD_PER_DEG (M_PI/180.0)
#define DEG_PER_RAD (180.0/M_PI)

bool valid_gps(Glib::ustring& text)
{
    auto regex = Glib::Regex::create("\\A[+\\-]?(\\d+\\.\\d*|\\d*\\.\\d+|\\d+)([eE][+\\-]?\\d+)?, *[+\\-]?(\\d+\\.\\d*|\\d*\\.\\d+|\\d+)([eE][+\\-]?\\d+)?\\z");
    return regex->match(text);
}

bool valid_number(Glib::ustring& text)
{
    auto regex = Glib::Regex::create("\\A[+\\-]?(\\d+\\.\\d*|\\d*\\.\\d+|\\d+)([eE][+\\-]?\\d+)?\\z");
    return regex->match(text);
}

Glib::ustring gps_to_text(double latitude, double longitude)
{
    return Glib::ustring::compose("%1, %2", number_to_text(latitude*DEG_PER_RAD),
                                  number_to_text(longitude*DEG_PER_RAD));
}

Glib::ustring number_to_text(double value)
{
//    return Glib::ustring::format(std::setprecision(15), value);
    return Glib::ustring( Glib::Ascii::dtostr(value) );
}

void text_to_gps(Glib::ustring& text, double &latitude, double &longitude)
{
    auto regex = Glib::Regex::create(", *");
    auto strings = regex->split(text);
    double latitude_deg;
    double longitude_deg;
    text_to_number(strings[0], latitude_deg);
    text_to_number(strings[1], longitude_deg);
    latitude = latitude_deg * RAD_PER_DEG;
    longitude = longitude_deg * RAD_PER_DEG;
}

void text_to_number(Glib::ustring &text, double &value)
{
    value = atof(text.data());
}
