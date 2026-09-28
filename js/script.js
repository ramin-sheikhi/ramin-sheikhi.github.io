"use strict";

document.addEventListener("DOMContentLoaded", () => {
  const $ = (sel, ctx = document) => ctx.querySelector(sel);
  const $$ = (sel, ctx = document) => [...ctx.querySelectorAll(sel)];

  /* ---------- Theme toggle (remembered per browser) ---------- */
  const root = document.documentElement;
  const themeToggle = $("#themeToggle");

  const savedTheme = safeStorage("get", "theme");
  const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
  root.setAttribute("data-theme", savedTheme || (prefersDark ? "dark" : "light"));

  themeToggle.addEventListener("click", () => {
    const next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
    root.setAttribute("data-theme", next);
    safeStorage("set", "theme", next);
  });

  function safeStorage(action, key, value) {
    try {
      return action === "get" ? localStorage.getItem(key) : localStorage.setItem(key, value);
    } catch {
      return null;
    }
  }

  /* ---------- Mobile menu ---------- */
  const navToggle = $("#navToggle");
  const navMenu = $("#navMenu");

  const setMenu = (open) => {
    navMenu.classList.toggle("open", open);
    navToggle.classList.toggle("open", open);
    navToggle.setAttribute("aria-expanded", String(open));
  };

  navToggle.addEventListener("click", () => setMenu(!navMenu.classList.contains("open")));
  $$(".nav__link").forEach((link) => link.addEventListener("click", () => setMenu(false)));

  /* ---------- Header shadow, back-to-top, active link ---------- */
  const header = $("#header");
  const toTop = $("#toTop");
  const sections = $$("main section[id]");
  const navLinks = $$(".nav__link");

  const onScroll = () => {
    const y = window.scrollY;
    header.classList.toggle("scrolled", y > 20);
    toTop.classList.toggle("show", y > 500);

    let current = sections[0].id;
    sections.forEach((sec) => {
      if (y >= sec.offsetTop - window.innerHeight / 3) current = sec.id;
    });
    navLinks.forEach((link) =>
      link.classList.toggle("active", link.getAttribute("href") === `#${current}`)
    );
  };

  window.addEventListener("scroll", onScroll, { passive: true });
  onScroll();

  toTop.addEventListener("click", () => window.scrollTo({ top: 0, behavior: "smooth" }));

  /* ---------- Typing effect ---------- */
  const typedEl = $("#typed");
  const words = [
    "fire following earthquake.",
    "multi-hazard risk.",
    "debris-flow impacts.",
    "seismic resilience.",
  ];
  let wordIndex = 0;
  let charIndex = 0;
  let deleting = false;

  (function type() {
    const word = words[wordIndex];
    typedEl.textContent = word.slice(0, charIndex);

    if (!deleting && charIndex < word.length) {
      charIndex++;
      setTimeout(type, 90);
    } else if (!deleting) {
      deleting = true;
      setTimeout(type, 1600);
    } else if (charIndex > 0) {
      charIndex--;
      setTimeout(type, 45);
    } else {
      deleting = false;
      wordIndex = (wordIndex + 1) % words.length;
      setTimeout(type, 300);
    }
  })();

  /* ---------- Animated counters ---------- */
  const animateCounter = (el) => {
    const target = Number(el.dataset.target);
    const duration = 1600;
    const start = performance.now();

    const step = (now) => {
      const progress = Math.min((now - start) / duration, 1);
      const eased = 1 - Math.pow(1 - progress, 3);
      el.textContent = Math.round(target * eased) + (progress === 1 ? el.dataset.suffix || "" : "");
      if (progress < 1) requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  };

  /* ---------- Scroll reveal, language bars, counters ---------- */
  const observer = new IntersectionObserver(
    (entries, obs) => {
      entries.forEach((entry) => {
        if (!entry.isIntersecting) return;
        const el = entry.target;

        if (el.classList.contains("reveal")) el.classList.add("visible");
        if (el.classList.contains("lang__fill")) el.style.width = `${el.dataset.width}%`;
        if (el.classList.contains("counter")) animateCounter(el);

        obs.unobserve(el);
      });
    },
    { threshold: 0.15 }
  );

  $$(".reveal, .lang__fill, .counter").forEach((el) => observer.observe(el));

  /* ---------- Publication filter ---------- */
  const filters = $$(".filter");
  const pubs = $$(".pub");

  filters.forEach((btn) => {
    btn.addEventListener("click", () => {
      filters.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");

      const category = btn.dataset.filter;
      pubs.forEach((pub) => {
        const match = category === "all" || pub.dataset.category === category;
        pub.classList.toggle("hide", !match);
        if (match) pub.classList.add("visible");
      });
    });
  });

  /* ---------- Contact form validation ---------- */
  const form = $("#contactForm");
  const success = $("#formSuccess");
  const emailPattern = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;

  const validators = {
    name: (v) => (v.length >= 2 ? "" : "Please enter your name."),
    email: (v) => (emailPattern.test(v) ? "" : "Please enter a valid email."),
    message: (v) => (v.length >= 10 ? "" : "Message should be at least 10 characters."),
  };

  const validateField = (field) => {
    const error = validators[field.name](field.value.trim());
    const group = field.closest(".form__group");
    group.classList.toggle("error", Boolean(error));
    $(".form__error", group).textContent = error;
    return !error;
  };

  $$("input, textarea", form).forEach((field) =>
    field.addEventListener("input", () => validateField(field))
  );

  form.addEventListener("submit", (e) => {
    e.preventDefault();
    const fields = $$("input, textarea", form);
    const allValid = fields.map(validateField).every(Boolean);
    if (!allValid) {
      success.textContent = "";
      return;
    }

    const button = $("button[type=submit]", form);
    button.disabled = true;
    button.textContent = "Sending…";

    // Simulated send — replace with a real API call (e.g. fetch) when you have a backend.
    setTimeout(() => {
      form.reset();
      button.disabled = false;
      button.textContent = "Send message";
      success.textContent = "Thanks! Your message has been sent.";
    }, 1200);
  });

  /* ---------- Footer year ---------- */
  $("#year").textContent = new Date().getFullYear();
});
