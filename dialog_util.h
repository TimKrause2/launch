#include <gtkmm.h>

bool valid_gps(Glib::ustring& text);
bool valid_number(Glib::ustring& text);

Glib::ustring gps_to_text(double latitude, double longitude);
Glib::ustring number_to_text(double value);
void text_to_gps(Glib::ustring& text, double &latitude, double &longitude);
void text_to_number(Glib::ustring& text, double &value);

Glib::RefPtr<Gtk::CssProvider> new_entry_provider(void);
void entry_set_font(Gtk::Widget &widget,
                    Glib::RefPtr<Gtk::CssProvider> provider);
