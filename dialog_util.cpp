#include "dialog_util.h"
#include <ios>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <cmath>


#define RAD_PER_DEG (M_PI/180.0)
#define DEG_PER_RAD (180.0/M_PI)
Glib::ustring gpsAlertMsg("Unrecognized GPS coordinate: ");
Glib::ustring numberAlertMsg("Unrecognized number: ");

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

void gpsAlertDialog(Gtk::Window &parent, Glib::ustring text)
{
    auto dialog = Gtk::AlertDialog::create(gpsAlertMsg + text);
    dialog->set_modal(true);
    dialog->show(parent);
}

void numberAlertDialog(Gtk::Window &parent, Glib::ustring text)
{
    auto dialog = Gtk::AlertDialog::create(numberAlertMsg + text);
    dialog->set_modal(true);
    dialog->show(parent);
}




Glib::RefPtr<Gtk::CssProvider> new_entry_provider(void)
{
    Glib::RefPtr<Gtk::CssProvider> provider = Gtk::CssProvider::create();
    provider->load_from_data(
                "#Entry {\n"
                "font-family:Monospace;\n"
                "font-size:18pt;\n"
                "}\n");
    return provider;
}

void entry_set_font(Gtk::Widget &widget,
                    Glib::RefPtr<Gtk::CssProvider> provider)
{
    widget.set_name("Entry");
    widget.get_style_context()->add_provider(provider, 1);
}
