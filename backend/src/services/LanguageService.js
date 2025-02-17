const translations = require('../data/translations.json');

class LanguageService {
    constructor() {
        this.defaultLanguage = 'he';
        this.supportedLanguages = ['he', 'en'];
    }

    getTranslation(key, language = this.defaultLanguage) {
        if (!this.supportedLanguages.includes(language)) {
            language = this.defaultLanguage;
        }

        const keys = key.split('.');
        let translation = translations[language];

        for (const k of keys) {
            if (!translation || !translation[k]) {
                return key; // Return the key if translation is not found
            }
            translation = translation[k];
        }

        return translation;
    }

    translateCourse(course, language = this.defaultLanguage) {
        if (!course) return null;

        return {
            ...course,
            name: course.name ? course.name[language] || course.name[this.defaultLanguage] : course.id,
            difficulty: this.getTranslation(`levels.${course.difficulty}`, language),
            category: this.getTranslation(`categories.${course.category}`, language)
        };
    }

    translateCourseList(courses, language = this.defaultLanguage) {
        return courses.map(course => this.translateCourse(course, language));
    }

    getSupportedLanguages() {
        return this.supportedLanguages;
    }

    getDefaultLanguage() {
        return this.defaultLanguage;
    }
}

module.exports = new LanguageService(); 