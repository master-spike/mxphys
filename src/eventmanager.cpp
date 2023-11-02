
#include "mxphys/event/eventmanager.h"

namespace mxphys::event {

void eventmanager::add_event(event_type etype, std::unique_ptr<event_data>&& edata) {
    auto it = m_eventdata.find(etype);
    if (it == m_eventdata.end()) {
        it = m_eventdata.emplace(etype, std::vector<std::unique_ptr<event_data>>()).first;
    }
    it->second.emplace_back(std::move(edata));
}
void eventmanager::emit_all() {
    for(auto& [etype, listeners] : m_registered_listeners) {
        auto data_it = m_eventdata.find(etype);
        if (data_it == m_eventdata.end()) continue;
        auto const& data_vec = data_it->second;
        for (auto& listener : listeners) {
            for (auto const& edata : data_vec) {
                listener(edata);
            }
        }
    }
    m_eventdata.clear();
}
void eventmanager::register_listener(event_type etype, const event_listener& listener) {
    auto it = m_registered_listeners.find(etype);
    if (it == m_registered_listeners.end()) {
        it = m_registered_listeners.emplace(etype, std::vector<event_listener>()).first;
    }
    it->second.emplace_back(listener);
}

}